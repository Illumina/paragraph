
PROTOBUF_GENERATE_CPP(PROTO_SRCS PROTO_HDRS graph.proto read.proto variant.proto)

add_library(proto STATIC ${PROTO_SRCS} ${PROTO_HDRS})
set (PROTO_LIBRARY proto PARENT_SCOPE)
