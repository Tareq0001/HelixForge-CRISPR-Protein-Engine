"""
HelixForge 3D - Smith-Waterman & Needleman-Wunsch Sequence Alignment Engine
Author: Tareq Aboushi
"""

from typing import Tuple, List


def smith_waterman_align(seq_a: str, seq_b: str, match: int = 3, mismatch: int = -2, gap: int = -2) -> Tuple[str, str, int, List[List[int]]]:
    m, n = len(seq_a), len(seq_b)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    max_score = 0
    max_pos = (0, 0)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s = match if seq_a[i - 1] == seq_b[j - 1] else mismatch
            val = max(0, dp[i - 1][j - 1] + s, dp[i - 1][j] + gap, dp[i][j - 1] + gap)
            dp[i][j] = val
            if val > max_score:
                max_score = val
                max_pos = (i, j)

    align_a, align_b = [], []
    i, j = max_pos
    while i > 0 and j > 0 and dp[i][j] > 0:
        s = match if seq_a[i - 1] == seq_b[j - 1] else mismatch
        if dp[i][j] == dp[i - 1][j - 1] + s:
            align_a.append(seq_a[i - 1])
            align_b.append(seq_b[j - 1])
            i -= 1
            j -= 1
        elif dp[i][j] == dp[i - 1][j] + gap:
            align_a.append(seq_a[i - 1])
            align_b.append("-")
            i -= 1
        else:
            align_a.append("-")
            align_b.append(seq_b[j - 1])
            j -= 1

    return "".join(reversed(align_a)), "".join(reversed(align_b)), max_score, dp
