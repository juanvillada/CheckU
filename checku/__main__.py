"""
Entry point for running CheckU as a module.

This allows running CheckU with:
python -m checku [commands]
"""

from .checku import app

if __name__ == "__main__":
    app()
