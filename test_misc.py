import re

input_str = 'ins(21)(p11.2p12p12). duplicated-insertion of Chr21: 6,239,312-5,700,622 (p12 - p12) into Chr21: 8,825,990 (p11.2)'
pattern = r'Chr(\d+): (\d{1,3}(?:,\d{3})*) \(.*?\)'
matches_itr = re.finditer(pattern, input_str)
for match in matches_itr:
    replacement_str = input_str[match.start(): match.end()]
    start_pos = int(match.group(2).replace(',', ''))
    print('chr' + match.group(1))
    print(start_pos)
