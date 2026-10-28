"""Shared paths for the private journal checks; artifacts stay outside sources."""
import os
from pathlib import Path
import tempfile

ROOT = Path(__file__).resolve().parents[4]
BUILD = ROOT/'siqs'
RESULTS = Path(os.environ.get('SIQS_JOURNAL_TEST_DIR') or
               tempfile.mkdtemp(prefix='siqs-journal-check-')).resolve()
RESULTS.mkdir(parents=True, exist_ok=True)
EXE = Path(os.environ.get('SIQS_JOURNAL_TEST_EXECUTABLE', str(BUILD/'msiqs')))
TOOL = BUILD/'tools/journal-info.py'
