from collections import defaultdict

d = defaultdict(list)
d[1].append('hi')
d[1].append('world')
d[2].append('yo')

print(d)