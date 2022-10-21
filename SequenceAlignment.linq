<Query Kind="Program" />

void Main()
{
	int mismatchPenalty = 5, gapPenalty = 4;
	var seq1 = "ACACATGCATCATGACTATGCATGCATGACTGACTGCATGCATGCATCCATCATGCATGCATCGATGCATGCATGACCACCTGTGTGACACATGCATGCGTGTGACATGCGAGACTCACTAGCGATGCATGCATGCATGCATGCATGC";
	var seq2 = "ATGATCATGCATGCATGCATCACACTGTGCATCAGAGAGAGCTCTCAGCAGACCACACACACGTGTGCAGAGAGCATGCATGCATGCATGCATGCATGGTAGCTGCATGCTATGAGCATGCAG";
	var nw = NeedlemanWunschSequenceAlignment(seq1, seq2, mismatchPenalty, gapPenalty);
	nw.Dump();
	VerifyScore(nw, mismatchPenalty, gapPenalty).Dump();
}

int VerifyScore((int, string[]) nwResult, int mismatchPenalty, int gapPenalty)
{
	var verifiedScore = 0;
	for (var i = 0; i < nwResult.Item2[0].Length; i++)
	{
		var e1 = nwResult.Item2[0][i];
		var e2 = nwResult.Item2[1][i];
		if (e1 != e2)
		{
			if (e1 == '-' || e2 == '-')
			{
				verifiedScore += gapPenalty;
			}
			else
			{
				verifiedScore += mismatchPenalty;
			}
		}
	}
	
	return verifiedScore;
}

(int, string[]) NeedlemanWunschSequenceAlignment(string seq1, string seq2, int mismatchPenalty, int gapPenalty)
{
	var cost = SequenceAlignmentPenalties(seq1, seq2, mismatchPenalty, gapPenalty, out int[][] subProblems);
	
	var s1Pos = seq1.Length;
	var s2Pos = seq2.Length;
	var al1 = new StringBuilder();
	var al2 = new StringBuilder();
	while (s1Pos > 0 && s2Pos > 0)
	{
		var diag = subProblems[s1Pos - 1][s2Pos - 1];
		var top = subProblems[s1Pos - 1][s2Pos];
		var left = subProblems[s1Pos][s2Pos - 1];
		var min = new[] { diag, top, left }.Min();
		if (seq1[s1Pos - 1] != seq2[s2Pos - 1])
		{
			if (left == min)
			{
				al1.Append('-');
				al2.Append(seq2[--s2Pos]);
			}
			else if (top == min)
			{
				al1.Append(seq1[--s1Pos]);
				al2.Append('-');
			}
			else
			{
				al1.Append(seq1[--s1Pos]);
				al2.Append(seq2[--s2Pos]);
			}
		}
		else
		{
			al1.Append(seq1[--s1Pos]);
			al2.Append(seq2[--s2Pos]);
		}
	}
	while (s1Pos-- > 0)
	{
		al1.Append(seq1[s1Pos]);
		al2.Append('-');
	}
	while (s2Pos-- > 0)
	{
		al2.Append(seq2[s2Pos]);
		al1.Append('-');
	}
	return (cost, new[] { new string(al1.ToString().Reverse().ToArray()), new string(al2.ToString().Reverse().ToArray()) });
}

int SequenceAlignmentPenalties(string seq1, string seq2, int mismatchPenalty, int gapPenalty, out int[][] subProblems)
{
	subProblems = new int[seq1.Length + 1][];
	for (var i = 0; i <= seq1.Length; i++)
	{
		subProblems[i] = new int[seq2.Length + 1];
		subProblems[i][0] = i * gapPenalty;
	}
	for (var j = 0; j <= seq2.Length; j++)
	{
		subProblems[0][j] = j * gapPenalty;
	}
	for (var i = 1; i <= seq1.Length; i++)
	{
		for (var j = 1; j <= seq2.Length; j++)
		{
			var diag = subProblems[i - 1][j - 1];
			var top = subProblems[i - 1][j] + gapPenalty;
			var left = subProblems[i][j - 1] + gapPenalty;
			var isMatch = seq1[i - 1] == seq2[j - 1];
			var min = new[] { diag + (isMatch ? 0 : mismatchPenalty), top, left }.Min();
			subProblems[i][j] = min;
		}
	}
	return subProblems[seq1.Length][seq2.Length];
}