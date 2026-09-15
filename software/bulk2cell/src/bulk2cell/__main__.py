"""Enable `python -m bulk2cell` using the installed command-line entry point."""
from .cli import main
raise SystemExit(main())
