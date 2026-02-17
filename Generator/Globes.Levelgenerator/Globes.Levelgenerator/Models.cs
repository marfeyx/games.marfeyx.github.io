namespace LevelGen;

public sealed class LevelObject
{
    public int Version { get; set; } = 1;
    public int N { get; set; } = 10;
    public int K { get; set; } = 2;

    public int Level { get; set; }

    // Estimated difficulty (from brute-force nodes)
    public int Difficulty { get; set; }

    public bool Unique { get; set; }
    public bool HumanSolvable { get; set; }
    public bool EasyStart { get; set; }

    public required int[][] Regions { get; set; }
    public required int[][] Solution { get; set; }
}

public sealed class SolvabilityInfo
{
    public required bool Solved { get; set; }
    public required bool EasyStart { get; set; }
    public required int ForcedFirst { get; set; }
    public required int Solutions { get; set; }
    public required int Nodes { get; set; }
    public required string FailReason { get; set; }
}

