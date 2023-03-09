#!/usr/bin/env python3
from .cli import cli 


def main() -> None:
    """
    Driver code for COMET.
    """
    args = cli()
    print(args)
    exit(0)


if __name__ == "__main__":
    main()
