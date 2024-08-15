import sys, json
slice_data = json.load(sys.stdin)['SliceTiming']
for s in slice_data: print s
