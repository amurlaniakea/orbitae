# Copyright (C) Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
from orbitae.atoms import BIO, BIO_EN, IONS, ION_BIO, ION_BIO_EN, NAMES, NAMES_EN, atom_record


def test_element_names_have_an_english_counterpart():
    assert len(NAMES) == len(NAMES_EN) == 36 and all(NAMES_EN)
    assert NAMES_EN[25] == "Iron" and NAMES_EN[10] == "Sodium"


def test_every_biomedical_text_is_translated():
    assert set(BIO) == set(BIO_EN) and set(ION_BIO) == set(ION_BIO_EN)
    for es, en in list(zip(BIO.values(), BIO_EN.values())) + list(zip(ION_BIO.values(), ION_BIO_EN.values())):
        assert en and en != es


def test_records_carry_both_languages():
    for Z, q in [(26, 0), (26, 2), (17, -1)] + IONS:
        rec = atom_record(Z, q)
        assert rec["name_en"] == NAMES_EN[Z - 1]
        assert bool(rec["bio"]) == bool(rec["bio_en"])
