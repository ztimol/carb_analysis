import os
import numpy as np


def is_whole_string_comment(s):
    """ function to check if a line
         starts with some character.
         Here # for comment
    """
    # return true if a line starts with #
    return s.lstrip().startswith("#")


def is_empty_string(string):
    if string.strip() in ["\n", "\r\n", ""]:
        return True
    return False


def clean_string(string):
    string = remove_string_comments(string)
    # return string.lstrip().rstrip().replace("\r", "").replace("\n", "")
    return string.strip()


def remove_string_comments(string):
    if string.find("#") != -1:
        return string[: string.find("#")]
    return string


def calc_exp_polynomial(x, a, b):
    return a * np.exp(b * x)


def calc_first_degree_polynomial(x, a, b):
    return a * x + b


def calc_second_degree_polynomial(x, a, b, c):
    return a * x ** 2 + b * x + c


def calc_third_order_polynomial(x, a, b, c, d):
    return a * x ** 3 + b * x ** 2 + c * x + d
