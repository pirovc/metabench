def join_params(params):
    return "_".join(params)


def params2str(arg, val):
    return str(arg) + "=" + str(val)


def str2params(s):
    return " ".join(s.replace("=", " ").split("_"))


def json_wildcards(wdict):
    json_str = ""
    for p, v in wdict.items():
        json_str += '\\"' + p + '\\": \\"' + v + '\\",\n'
    return json_str
