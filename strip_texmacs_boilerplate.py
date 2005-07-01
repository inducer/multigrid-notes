import sys
import re
import itertools

lines = file(sys.argv[1], "r").readlines()
lines = list(itertools.dropwhile(lambda x: x.find("\section{") == -1, lines))
lines = filter(lambda x: x.find("\end{doc") == -1, lines)
lines = [re.sub(r"\\epsfig\{file=([-a-z0-9]*)\.fig\}", r"\\input{\1.pstex_t}", line) for line in lines]

outfile = file(sys.argv[2], "w")
for line in lines:
    outfile.write(line)
