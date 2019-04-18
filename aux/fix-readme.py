#!/usr/bin/env python3
import re
from subprocess import Popen, PIPE

proc = Popen(['mhpl8r', '--help'], stdout=PIPE, stderr=PIPE, universal_newlines=True)
out, err = proc.communicate()

with open('README.md', 'r') as fh:
    text = fh.read().strip()
text = re.sub(
    r'..-- replaceme:start --..*..-- replaceme.end --.',
    '<!-- replaceme:start -->\n```\n' + out +'\n```\n<!-- replaceme.end -->',
    text,
    flags=re.DOTALL
)
with open('README.md', 'w') as fh:
    print(text, file=fh)
