# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Camille Clouard"
__copyright__ = "Copyright 2024, Camille Clouard"
__email__ = "camille.clouard@scilifelab.uu.se"
__license__ = "GPL-3"


def test_dummy():
    from dummy import dummy
    assert dummy() == 1
