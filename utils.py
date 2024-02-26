def sign(x):
    if x < 0:
        return -1
    elif x == 0:
        return 0
    else:
        return 1


def reverse_dict(input_dict):
    """
    requires mapping to be bijective
    :param input_dict:
    :return:
    """
    output_dict = {}
    for position, name in input_dict.items():
        output_dict[name] = position
    return output_dict
