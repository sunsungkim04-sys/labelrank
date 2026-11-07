import sys

from .cli import main

try:
    sys.exit(main())
except BrokenPipeError:   # output piped into head and the like
    sys.stderr.close()
    sys.exit(141)
