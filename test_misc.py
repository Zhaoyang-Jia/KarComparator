# import re
#
# input_str = 'ins(21)(p11.2p12p12). duplicated-insertion of Chr21: 6,239,312-5,700,622 (p12 - p12) into Chr21: 8,825,990 (p11.2)'
# pattern = r'Chr(\d+): (\d{1,3}(?:,\d{3})*) \(.*?\)'
# matches_itr = re.finditer(pattern, input_str)
# for match in matches_itr:
#     replacement_str = input_str[match.start(): match.end()]
#     start_pos = int(match.group(2).replace(',', ''))
#     print('chr' + match.group(1))
#     print(start_pos)

# import pandas as pd
#
# # # Sample DataFrame
# # data = {'file_name': ['file1', 'file1', 'file2', 'file2', 'file3'],
# #         'karsim_CN': [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0], [10.0, 11.0, 12.0], [13.0, 14.0, 15.0]]}
# # df = pd.DataFrame(data)
# #
# # def vector_addition(lists):
# #     result = [sum(x) for x in zip(*lists)]
# #     return result
# #
# # sums_by_file = df.groupby('file_name')['karsim_CN'].apply(vector_addition)
# #
# # # Convert the resulting Series to a dictionary
# # result_dict = sums_by_file.to_dict()
# #
# # # Display the resulting dictionary
# # print(result_dict)
# import numpy as np
# def sum_of_absolute_differences(array1, array2):
#     # Compute the absolute differences between the arrays
#     absolute_differences = np.abs(array1 - array2)
#
#     # Sum the absolute differences
#     sum_absolute_differences = np.sum(absolute_differences)
#
#     return sum_absolute_differences
#
#
# # Example arrays
# array1 = np.array([1, 5, 3])
# array2 = np.array([4, 2, 6])
#
# # Calculate the sum of the absolute value of the differences
# result = sum_of_absolute_differences(array1, array2)
#
# print("Sum of absolute differences:", result)

# import base64
# img = "/Users/zhaoyangjia/PyCharm_Repos/KarComparator/latex_reports/paul_dremsek_plots_new/3_imagecluster0_rotated.png"
#
# def image_to_base64(image_path):
#     try:
#         with open(image_path, "rb") as img_file:
#             return base64.b64encode(img_file.read()).decode('utf-8')
#     except FileNotFoundError:
#         print(f"Error: File {image_path} not found.")
#         return ""
#
# print(image_to_base64(img))

# x = ['a', 'b', 'c']
# y = ['d', 'e', 'f']
# for z in x + y:
#     print(z)

# import re
#
# # Example string
# text = "This is an example (with some text inside parentheses) and (another example)."
#
# # Regular expression to match everything inside parentheses
# pattern = r'\(.*?\)'
#
# # Find all matches in the text
# matches = re.findall(pattern, text)
#
# print(matches)

# x = [1, 2, 3]
# y = {1, 5, 6}
#
# z = y.intersection(x)
# print(len(z))

import re
# match = re.match(r'(\d+)(\D*)', '0a')
# print(match.group(1))
# print(match.group(2))
# print(ord(match.group(2)) - 96)
# print()

def split_int_char(input_string):
    # Regular expression to match an integer followed by 0 to 3 characters
    match = re.match(r'(\d+)(\D*)', input_string)
    if match:
        int_part = match.group(1)
        char_part = match.group(2)
        if len(match.group(2)) == 0:
            char_value = 0.0
        else:
            char_value = 0.0
            multiplier = 0.01
            for character in match.group(2):
                char_value += multiplier * (ord(character) - 96)
                multiplier *= 0.01
        return int_part, char_part, int(match.group(1)) + char_value
    return None, None, None  # In case there's no match

# Example usage
input_string1 = "123a"
input_string2 = "456ab"
input_string3 = "789abc"

int_part1, char_part1, v1 = split_int_char(input_string1)
int_part2, char_part2, v2 = split_int_char(input_string2)
int_part3, char_part3, v3 = split_int_char(input_string3)

print(f"Input: {input_string1} -> Integer part: {int_part1}, Character part: {char_part1}, value: {v1}")
print(f"Input: {input_string2} -> Integer part: {int_part2}, Character part: {char_part2}, value: {v2}")
print(f"Input: {input_string3} -> Integer part: {int_part3}, Character part: {char_part3}, value: {v3}")
