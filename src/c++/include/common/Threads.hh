// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2017 Illumina, Inc.
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/**
 * \brief Threading helpers
 *
 * \file Threads.hh
 * \author Roman Petrovski
 * \email rpetrovski@illumina.com
 *
 */

#pragma once

#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>

#include "Error.hh"

namespace common
{

namespace detail
{

    /**
     * \brief Base non-template class for ScopeEndCall
     */
    class ScopeEndCallBase
    {
    public:
        virtual ~ScopeEndCallBase() = default;

        explicit operator bool() const { return false; }
    };

    /**
     * \brief Holds a functor which is called at the destruction time
     */
    template <typename FuncT> class ScopeEndCall : public ScopeEndCallBase
    {
        FuncT f_;

    public:
        explicit ScopeEndCall(FuncT f)
            : f_(f)
        {
        }

        ~ScopeEndCall() override { f_(std::uncaught_exception()); }
    };

    /**
     * \brief Helper to create ScopeEndCallHolder
     */
    template <typename FuncT> const detail::ScopeEndCall<FuncT> makeScopeEndCallHolder(FuncT f)
    {
        return detail::ScopeEndCall<FuncT>(f);
    }
}

/**
 * \brief ensures f is called during the stack unwind of the scope following the macro
 */
#define ASYNC_BLOCK_WITH_CLEANUP(f)                                                                                    \
    if (const common::detail::ScopeEndCallBase& b = common::detail::makeScopeEndCallHolder(f))                         \
    {                                                                                                                  \
        (void)b;                                                                                                       \
    }                                                                                                                  \
    else

/**
 * \brief Inversion of the std::unique_lock and such
 */
template <typename Lock> class unlock_guard
{
private:
    Lock& l;

public:
    unlock_guard(unlock_guard&) = delete;
    unlock_guard& operator=(unlock_guard&) = delete;
    explicit unlock_guard(Lock& m_)
        : l(m_)
    {
        l.unlock();
    }

    ~unlock_guard() { l.lock(); }
};

template <bool crashOnExceptions> class BasicThreadPool
{
    struct Executor
    {
        Executor(const std::size_t maxThreads, const unsigned request)
            : maxThreads_(maxThreads)
            , request_(request)
        {
        }
        virtual void execute() = 0;
        virtual ~Executor() = default;

        const std::size_t maxThreads_;
        const int request_;
        Executor* next_ = 0;
        std::size_t threadsIn_ = 0;
        bool complete_ = false;
        friend std::ostream& operator<<(std::ostream& os, const Executor& e)
        {
            return os << "Executor(" << e.maxThreads_ << "mt," << e.request_ << "r," << e.threadsIn_ << "ti,"
                      << e.complete_ << "c)";
        }
    } * head_;
    static std::mutex mutex_;
    std::condition_variable stateChangedCondition_;

    // true when the whole thing goes down
    bool terminateRequested_;

    // number of threads that have entered ready state
    std::size_t threadsReady_ = 0;

    // constantly incrementing number to make sure each thread processes one master call only once
    static int CURRENT_REQUEST_;
    static int __thread THREAD_LAST_PROCESSED_REQUEST_;

    std::exception_ptr firstThreadException_;

    typedef std::vector<std::thread> ThreadVector;
    typedef ThreadVector::size_type size_type;
    ThreadVector threads_{};

public:
    /**
     * \return number of threads in th pool + 1. This is because the thread calling execute()
     *         is counted as a worker.
     */
    std::size_t size() const { return threads_.size() + 1; }

    /**
     * \brief clears and repopulates ThreadPool. Use at your own risk. The only reasonable situation you'd might
     *        want to do this is in the sequence of unit tests that require different parallelization.
     *        also resets terminateRequested_ state.
     */
    void reset(std::size_t newSize)
    {
        clear();
        std::unique_lock<std::mutex> lock(mutex_);
        assert(0 < newSize); //, "Inadequate pool size";
        // thread calling the execute will be one of the workers
        while (--newSize)
        {
            threads_.push_back(std::thread(&BasicThreadPool::threadFunc, this));
        }
    }

    /**
     * \brief Constructs a vector of size threads. All memory allocations that are required happen
     *        at this point.
     *
     *  @param size         number of threads to produce
     *  @param numaNode     NUMA node or one of the numa::default* constants
     */
    explicit BasicThreadPool(size_type size)
        : head_(0)
        , terminateRequested_(false)
    {
        reset(size);
    }

    /**
     * \brief Tells all threads to terminate and releases them.
     */
    ~BasicThreadPool() { clear(); }

    /**
     * \brief Executes func on requested number of threads.
     *
     * \threads number of threads to use. This must be less or equal than size()
     * \param func Unary function to execute.
     */
    template <typename F> void execute(F func, const unsigned threads)
    {
        std::unique_lock<std::mutex> lock(mutex_);
        if (firstThreadException_)
        {
            LOG()->warn("WARNING: execute called when an exception is pending ");
            std::rethrow_exception(firstThreadException_);
        }

        assert(threads <= size()); //, "Request must not exceed the amount of threads available");
        struct FuncExecutor : public Executor
        {
            F& func_;
            FuncExecutor(F& func, const unsigned threads)
                : Executor(threads, ++CURRENT_REQUEST_)
                , func_(func)
            {
            }
            virtual void execute() { func_(); }
        } executor(func, threads);

        LOG()->trace("created {}", executor);

        enque(&executor);
        LOG()->trace("enqued {}", executor);
        stateChangedCondition_.notify_all();

        while (!terminateRequested_ && (!executor.complete_ || executor.threadsIn_))
        {
            executeAllPending(lock);
            if (!terminateRequested_ && (!executor.complete_ || executor.threadsIn_))
            {
                stateChangedCondition_.wait(lock);
            }
        }

        unque(&executor);
        LOG()->trace("unqued {}", executor);

        if (firstThreadException_)
        {
            LOG()->warn("WARNING: rethrowing a thread exception ");
            std::rethrow_exception(firstThreadException_);
        }
    }

    /**
     * \brief Executes func on requested number of size() threads.
     **/
    template <typename F> void execute(F func) { execute(func, static_cast<const unsigned int>(size())); }

private:
    void clear()
    {
        {
            std::unique_lock<std::mutex> lock(mutex_);
            // if clear is called immediately after construction, some threads might not be ready to be joined
            waitAllReady(lock);

            terminateRequested_ = true;
            stateChangedCondition_.notify_all();
        }
        std::for_each(threads_.begin(), threads_.end(), [](std::thread& t) { t.join(); });
        threads_.clear();
        terminateRequested_ = false;
        threadsReady_ = 0;
    }

    void unque(Executor* executor)
    {
        assert(head_); //, "Expected executor to be set");
        Executor* prev = 0;
        Executor* e = head_;
        for (; executor != e; prev = e, e = e->next_)
        {
            ;
        }

        assert(executor == e);
        if (!prev)
        {
            head_ = executor->next_;
        }
        else
        {
            prev->next_ = executor->next_;
        }
    }

    void enque(Executor* executor)
    {
        executor->next_ = head_;
        head_ = executor;
    }

    /**
     * \brief find next suitable executor for the current thread
     */
    Executor* findNext()
    {
        Executor* ret = 0;
        for (Executor* e = head_; 0 != e; e = e->next_)
        {
            // LOG()->trace("checking {}", *e);

            if (!e->complete_ && e->threadsIn_ < e->maxThreads_ && THREAD_LAST_PROCESSED_REQUEST_ < e->request_
                && (!ret || ret->request_ > e->request_))
            {
                ret = e;
                // LOG()->trace("found {}", *e);
            }
        }
        if (!ret)
        {
            // LOG()->trace("found nothing");
        }
        return ret;
    }

    /**
     * \brief executes and allows the exception to escape
     */
    void unsafeExecute(std::unique_lock<std::mutex>& lock, Executor* executor)
    {
        assert(!executor->complete_);
        assert(executor->threadsIn_ < executor->maxThreads_);

        ++executor->threadsIn_;
        THREAD_LAST_PROCESSED_REQUEST_ = executor->request_;
        try
        {
            unlock_guard<std::unique_lock<std::mutex>> unlock(lock);
            executor->execute();
        }
        catch (...)
        {
            executor->complete_ = true;
            --executor->threadsIn_;
            stateChangedCondition_.notify_all();
            throw;
        }
        executor->complete_ = true;
        --executor->threadsIn_;
        stateChangedCondition_.notify_all();
    }

    /**
     * \brief executes and stores the exception information so that it can be rethrown on the main thread
     */
    void safeExecute(std::unique_lock<std::mutex>& lock, Executor* executor)
    {
        try
        {
            unsafeExecute(lock, executor);
        }
        catch (...)
        {
            if (!firstThreadException_)
            {
                firstThreadException_ = std::current_exception();
                LOG()->critical("ERROR: This thread caught an exception first");
            }
            else
            {
                LOG()->critical("ERROR: This thread also caught an exception");
            }
        }
    }

    void threadFunc()
    {
        std::unique_lock<std::mutex> lock(mutex_);
        ++threadsReady_;
        stateChangedCondition_.notify_all();
        while (!terminateRequested_)
        {
            executeAllPending(lock);
            if (!terminateRequested_)
            {
                stateChangedCondition_.wait(lock);
            }
        }
    }

    void waitAllReady(std::unique_lock<std::mutex>& lock)
    {
        while (size() != (threadsReady_ + 1))
        {
            stateChangedCondition_.wait(lock);
        }
    }

    void executeAllPending(std::unique_lock<std::mutex>& lock)
    {
        for (Executor* executor = findNext(); !terminateRequested_ && 0 != executor; executor = findNext())
        {
            if (crashOnExceptions)
            {
                unsafeExecute(lock, executor);
            }
            else
            {
                safeExecute(lock, executor);
            }
        }
    }
};
template <bool crashOnExceptions> std::mutex BasicThreadPool<crashOnExceptions>::mutex_;
template <bool crashOnExceptions> int BasicThreadPool<crashOnExceptions>::CURRENT_REQUEST_ = 0;
template <bool crashOnExceptions> int __thread BasicThreadPool<crashOnExceptions>::THREAD_LAST_PROCESSED_REQUEST_ = 0;

typedef BasicThreadPool<false> SafeThreadPool;
typedef BasicThreadPool<true> UnsafeThreadPool;
typedef SafeThreadPool ThreadPool;

/**
 * \param threadsMax  maximum number of threads. Used in the singleton initialization.
 *                    will produce failure if changes on subsequent calls.
 *                    special value 0 will result in pool allocation with std::thread::hardware_concurrency(),
 *                    also will not lead to a failure on subsequent calls regardless of the size of allocated pool
 * \return reference to a thread pool singleton
 */
ThreadPool& CPU_THREADS(std::size_t threadsMax = 0);

} // namespacecommon
