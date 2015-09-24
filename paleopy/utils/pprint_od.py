def pprint_od(od):
    print("{")
    for key in od:
        print("{}:{}".format(key, od[key]))
    print("}")