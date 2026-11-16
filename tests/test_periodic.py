# Copyright (C) Pedro Sordo Martínez
# SPDX-License-Identifier: AGPL-3.0-or-later
import pytest

from orbitae.atoms import CATEGORY, atom_record, group_of, period_of


@pytest.mark.parametrize("Z,g,p", [(1, 1, 1), (2, 18, 1), (5, 13, 2), (10, 18, 2), (13, 13, 3), (17, 17, 3),
                                   (19, 1, 4), (26, 8, 4), (30, 12, 4), (36, 18, 4)])
def test_group_and_period(Z, g, p):
    assert (group_of(Z), period_of(Z)) == (g, p)


def test_periods_have_the_right_sizes():
    sizes = [sum(period_of(Z) == p for Z in range(1, 37)) for p in (1, 2, 3, 4)]
    assert sizes == [2, 8, 8, 18]
    assert sorted({group_of(Z) for Z in range(1, 37)}) == list(range(1, 19))


def test_categories_cover_all_elements():
    assert set(CATEGORY) == set(range(1, 37))
    assert set(CATEGORY.values()) <= {"nonmetal", "noble", "alkali", "alkaline", "metalloid", "halogen", "post", "transition"}
    assert (CATEGORY[26], CATEGORY[17], CATEGORY[11], CATEGORY[14], CATEGORY[18]) == ("transition", "halogen", "alkali", "metalloid", "noble")
    assert all(CATEGORY[Z] == "noble" for Z in (2, 10, 18, 36))


@pytest.mark.parametrize("Z,q,sub,n_out", [(26, 0, "4s", 2), (17, 0, "3p", 7), (11, 0, "3s", 1), (2, 0, "1s", 2),
                                            (29, 0, "4s", 1), (26, 2, "3d", 14), (17, -1, "3p", 8)])
def test_outer_shell_data(Z, q, sub, n_out):
    rec = atom_record(Z, q)
    assert (rec["outer_subshell"], rec["outer_electrons"]) == (sub, n_out)
    assert rec["group"] == group_of(Z) and rec["period"] == period_of(Z) and rec["category"] == CATEGORY[Z]
