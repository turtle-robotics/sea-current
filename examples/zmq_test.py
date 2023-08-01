#!/usr/bin/env python3

import zmq
import json
import time

path_str = '''-0.5 0.5 -0.25 0.25
0 0
10 0
10 10'''

if __name__ == '__main__':
    ctx = zmq.Context()
    s = ctx.socket(zmq.REQ)
    # s.connect('tcp://localhost:8080')
    s.connect('ipc:///tmp/test')
    offset_time = time.perf_counter()
    s.send(bytes(path_str, 'utf-8'))
    idk = s.recv()
    new_time = time.perf_counter()
    print(new_time - offset_time)
    j = json.loads(idk)
    # print(j)
    with open("output.json", "w") as json_file:
        json.dump(j, json_file, indent=4)
