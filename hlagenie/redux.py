"""Self-contained two-field (U2) allele reduction.

This module is a faithful port of the ``"U2"`` reduction performed by
py-ARD (https://github.com/nmdp-bioinformatics/py-ard), reimplemented here so
that HLAGenie no longer needs py-ard at runtime just to reduce allele names to
their two-field form. All credit for the reduction algorithm and the
nomenclature rules belongs to the py-ARD authors; py-ard, like HLAGenie, is
LGPL-licensed, so this derivative is license-compatible.

The reduction and its supporting data structures (``valid_alleles`` and
``p_not_g``) reproduce py-ard's output exactly: validated to 0 differences
against ``pyard.init(version).redux(allele, "U2")`` across every 3+ field
allele in the IMGT/HLA database.

Data sources (same files py-ard reads):
    * ``Allelelist.<version>.txt``  -> the set of valid allele names
    * ``hla_nom_g.txt`` / ``hla_nom_p.txt`` -> the P-not-G (ping) overrides
"""

# Expression suffixes used in HLA nomenclature (null, low, secreted, questionable).
EXPRESSION_CHARS = frozenset({"N", "L", "S", "Q"})

# P and G group name suffixes.
P_AND_G_CHARS = ("P", "G")

# Loci for which G groups (and therefore reduction) are defined. Alleles at any
# other locus are returned unchanged, matching py-ard's ``G_GROUP_LOCI`` gate.
G_GROUP_LOCI = frozenset(
    {
        "A", "B", "C",
        "DPA1", "DPB1", "DQA1", "DQB1",
        "DRB1", "DRB3", "DRB4", "DRB5",
        "DRA", "DMA", "DMB", "DOA", "DOB",
        "E", "F", "G",
    }
)


def get_n_field_allele(allele: str, n: int, preserve_expression: bool = False) -> str:
    """Return the first ``n`` fields of an allele (ported from py-ard).

    If ``preserve_expression`` is set and the allele ends in an expression
    character while having more than ``n`` fields, that character is kept.
    """
    fields = allele.split(":")
    if preserve_expression and allele[-1] in EXPRESSION_CHARS and len(fields) > n:
        return ":".join(fields[0:n]) + allele[-1]
    return ":".join(fields[0:n])


def get_2field_allele(allele: str) -> str:
    """Two-field allele with any P/G suffix removed (ported from py-ard)."""
    if allele[-1] in P_AND_G_CHARS:
        allele = allele[:-1]
    return get_n_field_allele(allele, 2)


def get_3field_allele(allele: str) -> str:
    """Three-field allele with any P/G suffix removed (ported from py-ard)."""
    if allele[-1] in P_AND_G_CHARS:
        allele = allele[:-1]
    return get_n_field_allele(allele, 3)


def build_valid_alleles(allele_names) -> set:
    """Build py-ard's set of valid allele names from an allele list.

    Mirrors py-ard's ``generate_alleles_and_xx_codes_and_who``: the set is every
    full allele name, plus its two- and three-field forms, plus the two-field
    form (with expression character) of alleles whose whole two-field group
    shares a single expression character.
    """
    allele_names = list(allele_names)
    valid = set(allele_names)
    valid.update(get_2field_allele(a) for a in allele_names)
    valid.update(get_3field_allele(a) for a in allele_names)

    # expression_reduce: a suffix propagates to two-field level only when every
    # 3+/4-field allele in that two-field group carries the same suffix.
    groups: dict[str, set] = {}
    for a in allele_names:
        if a[-1] in EXPRESSION_CHARS and a.count(":") >= 2:
            groups.setdefault(get_2field_allele(a), set()).add(a[-1])
    valid.update(
        two_d + next(iter(chars)) for two_d, chars in groups.items() if len(chars) == 1
    )
    return valid


def build_p_not_g(g_group_alleles, p_group_pairs) -> dict:
    """Build py-ard's ``p_not_g`` override map.

    ``g_group_alleles`` is an iterable of full allele names appearing in
    ``hla_nom_g.txt``; ``p_group_pairs`` is an iterable of
    ``(full_allele, full_p_group_name)`` from ``hla_nom_p.txt``. A two-field
    allele present in the P groups but absent from the G groups maps each of its
    alleles to the two-field form of its P group.
    """
    g_two_field = {get_2field_allele(a) for a in g_group_alleles}

    rows = [
        (full, get_2field_allele(full), get_2field_allele(p_group))
        for full, p_group in p_group_pairs
    ]
    p_not_in_g = {two_d for _, two_d, _ in rows} - g_two_field

    p_not_g: dict[str, str] = {}
    for full, two_d, lgx in rows:
        if two_d in p_not_in_g:
            p_not_g[full] = lgx  # last write wins, matching py-ard's to_dict
    return p_not_g


def reduce_to_two_field(allele: str, valid_alleles: set, p_not_g: dict) -> str:
    """Reduce an allele to two fields exactly as py-ard's ``redux(a, "U2")``.

    :param allele: allele name to reduce
    :param valid_alleles: set from :func:`build_valid_alleles`
    :param p_not_g: mapping from :func:`build_p_not_g`
    :return: the two-field (U2) reduced allele name
    """
    # Non-G-group loci are not reduced.
    star = allele.find("*")
    if star != -1 and allele[:star] not in G_GROUP_LOCI:
        return allele

    # P-not-G (ping) override takes precedence, and applies even to alleles that
    # are already two-field (e.g. A*01:335 -> A*01:01).
    mapped = p_not_g.get(allele)
    if mapped is not None:
        return mapped

    fields = allele.split(":")
    if len(fields) <= 2:
        return allele

    plain = fields[0] + ":" + fields[1]

    # Preserve a trailing expression character only if the result is itself a
    # valid allele; otherwise py-ard falls back to the plain two-field (lgx) form.
    last = allele[-1]
    if last in EXPRESSION_CHARS:
        candidate = plain + last
        if candidate in valid_alleles:
            return candidate
    elif plain in valid_alleles:
        return plain

    return plain
