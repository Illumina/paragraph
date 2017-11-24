#! /usr/bin/env python
from __future__ import print_function

import os
import json
import jsonschema

base_dir = os.path.abspath(os.path.join(__file__, '..', '..', '..'))
schema_dir = os.path.join(base_dir, 'share', 'schema')
paragraph_dir = os.path.join(base_dir, 'data', 'paragraph')

input_schema = os.path.join(schema_dir, 'input_schema.json')
output_schema = os.path.join(schema_dir, 'output_schema.json')
input_to_validate = [os.path.join(paragraph_dir, 'simple', 'del-example-3.json'),
                     os.path.join(paragraph_dir, 'simple', 'del-example-4.json'),
                     os.path.join(paragraph_dir, 'simple', 'swap-example-1.json'),
                     os.path.join(paragraph_dir, 'simple', 'swap-example-2.json'),
                     os.path.join(paragraph_dir, 'simple', 'swap-example-5.json'),
                     os.path.join(paragraph_dir, 'long-del', 'chr4-21369091-21376907.json'),
                     os.path.join(paragraph_dir, 'pg-het-ins', 'pg-het-ins.json')]

# missing sample IDs in both outputs, invalid for testing
# output_to_validate = [os.path.join(paragraph_dir, 'long-del', 'chr4-21369091-21376907.paragraph.json'),
#                       os.path.join(paragraph_dir, 'pg-het-ins', 'pg-het-ins.result.json')]


with open(input_schema) as schema_f:
    in_schema = json.load(schema_f)

with open(output_schema) as schema_f:
    out_schema = json.load(schema_f)

jsonschema.Draft4Validator.check_schema(in_schema)
jsonschema.Draft4Validator.check_schema(out_schema)

# from https://github.com/Julian/jsonschema/issues/313
schema_dir = os.path.dirname(os.path.abspath(input_schema))

resolver = jsonschema.RefResolver(base_uri='file://' + schema_dir + '/', referrer=in_schema)

for spec_f in input_to_validate:
    print(spec_f)
    with open(spec_f) as spec:
        jsonschema.validate(json.load(spec), in_schema, resolver=resolver)

resolver = jsonschema.RefResolver(base_uri='file://' + schema_dir + '/', referrer=out_schema)

for spec_f in output_to_validate:
    print(spec_f)
    with open(spec_f) as spec:
        jsonschema.validate(json.load(spec), out_schema, resolver=resolver)

