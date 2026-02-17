using SkiaSharp;

namespace LevelGen;

public static class SolutionPngRenderer
{
    public static void SaveSolutionPng(LevelObject level, string filepath, int cellSize = 32, int padding = 16)
    {
        Directory.CreateDirectory(Path.GetDirectoryName(filepath)!);

        int n = level.N;
        int w = padding * 2 + n * cellSize;
        int h = padding * 2 + n * cellSize;

        using var surface = SKSurface.Create(new SKImageInfo(w, h, SKColorType.Rgba8888, SKAlphaType.Premul));
        var canvas = surface.Canvas;
        canvas.Clear(SKColors.White);

        // Grid lines
        using var gridPaint = new SKPaint
        {
            Color = new SKColor(220, 220, 220),
            StrokeWidth = 1,
            IsAntialias = true,
            Style = SKPaintStyle.Stroke
        };

        for (int i = 0; i <= n; i++)
        {
            int x = padding + i * cellSize;
            int y = padding + i * cellSize;

            canvas.DrawLine(x, padding, x, padding + n * cellSize, gridPaint);
            canvas.DrawLine(padding, y, padding + n * cellSize, y, gridPaint);
        }

        // Region borders
        using var regPaint = new SKPaint
        {
            Color = SKColors.Black,
            StrokeWidth = 2,
            IsAntialias = true,
            Style = SKPaintStyle.Stroke
        };

        // Outer border
        canvas.DrawRect(padding, padding, n * cellSize, n * cellSize, regPaint);

        var rgn = level.Regions;
        for (int r = 0; r < n; r++)
        {
            for (int c = 0; c < n; c++)
            {
                int x0 = padding + c * cellSize;
                int y0 = padding + r * cellSize;

                if (c + 1 < n && rgn[r][c] != rgn[r][c + 1])
                {
                    int x = x0 + cellSize;
                    canvas.DrawLine(x, y0, x, y0 + cellSize, regPaint);
                }

                if (r + 1 < n && rgn[r][c] != rgn[r + 1][c])
                {
                    int y = y0 + cellSize;
                    canvas.DrawLine(x0, y, x0 + cellSize, y, regPaint);
                }
            }
        }

        // Globes
        int inset = (int)(cellSize * 0.18);
        int d = cellSize - inset * 2;

        using var globeFill = new SKPaint
        {
            Color = new SKColor(11, 65, 205),
            IsAntialias = true,
            Style = SKPaintStyle.Fill
        };

        using var globeStroke = new SKPaint
        {
            Color = SKColors.Black,
            StrokeWidth = 1,
            IsAntialias = true,
            Style = SKPaintStyle.Stroke
        };

        var sol = level.Solution;
        for (int r = 0; r < n; r++)
        {
            for (int c = 0; c < n; c++)
            {
                if (sol[r][c] != 1) continue;

                float x = padding + c * cellSize + inset;
                float y = padding + r * cellSize + inset;

                canvas.DrawOval(x + d / 2f, y + d / 2f, d / 2f, d / 2f, globeFill);
                canvas.DrawOval(x + d / 2f, y + d / 2f, d / 2f, d / 2f, globeStroke);
            }
        }

        using var image = surface.Snapshot();
        using var data = image.Encode(SKEncodedImageFormat.Png, 100);

        using var fs = File.OpenWrite(filepath);
        data.SaveTo(fs);
    }
}
