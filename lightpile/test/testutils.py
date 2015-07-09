# -*- coding: utf-8 -*-

import sys
import os
from contextlib import contextmanager
from StringIO import StringIO
from tempfile import NamedTemporaryFile
import numpy as np

__author__ = "Richard Pfeifer"


@contextmanager
def captured_err():
    """
    Context with temporary redirected stderr. Returns this temporary file
    like object.
    """
    new_err = StringIO()
    old_err = sys.stderr
    try:
        sys.stderr = new_err
        yield sys.stderr
    finally:
        sys.stderr = old_err


@contextmanager
def captured_out():
    """
    Context with temporary redirected stdout and stderr. Returns tuple of
    temporary file like objects replacing (stdout, stderr).
    """
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout = new_out
        sys.stderr = new_err
        yield (sys.stdout, sys.stderr)
    finally:
        sys.stdout = old_out
        sys.stderr = old_err


@contextmanager
def temp_file_path(bytestring_content):
    """
    Supply a context with the filename of a temporary file with the given
    content. The file is closed and has to be opened.

    It is ensured that the file will be closed and removed.
    """
    f = NamedTemporaryFile(mode='w+b', prefix="lightpile", delete=False)
    try:
        f.write(bytestring_content)
        f.close()
        yield f.name
    finally:
        f.close()
        os.remove(f.name)


@contextmanager
def np_printoptions(*args, **kwargs):
    """
    Provide context with custom numpy printoptions.
    """
    original_printoptions = np.get_printoptions()
    np.set_printoptions(*args, **kwargs)
    yield
    np.set_printoptions(original_printoptions)
