import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";
import readline from "readline/promises";
import { stdin as input, stdout as output } from "process";

import { generateLevel } from "./generator-core.js";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const levelsDir = path.resolve(__dirname, "../levels");
if (!fs.existsSync(levelsDir)) fs.mkdirSync(levelsDir, { recursive: true });

function parseIntStrict(s) {
  const n = Number(String(s).trim());
  return Number.isInteger(n) ? n : NaN;
}

function validateRange(min, max) {
  if (!Number.isInteger(min) || !Number.isInteger(max)) {
    throw new Error("Difficulty values must be integers.");
  }
  if (min < 1 || max > 10 || min > max) {
    throw new Error("Difficulty range must be between 1 and 10, and min <= max.");
  }
}

function getNextLevelNumber(dirPath) {
  const files = fs
    .readdirSync(dirPath, { withFileTypes: true })
    .filter(d => d.isFile())
    .map(d => d.name);

  let maxLevel = 0;
  const reNew = /^Level\s+(\d+)\s+(\d+)\.json$/i;
  const reOld = /^lvl(\d+)\.json$/i;

  for (const name of files) {
    let m = name.match(reNew);
    if (m) {
      const n = Number(m[1]);
      if (Number.isFinite(n)) maxLevel = Math.max(maxLevel, n);
      continue;
    }
    m = name.match(reOld);
    if (m) {
      const n = Number(m[1]);
      if (Number.isFinite(n)) maxLevel = Math.max(maxLevel, n);
    }
  }

  return maxLevel + 1;
}

function sleep(ms) {
  return new Promise(res => setTimeout(res, ms));
}

function randomInt(min, max) {
  return Math.floor(Math.random() * (max - min + 1)) + min;
}

const rl = readline.createInterface({ input, output });

try {
  const amount = parseIntStrict(await rl.question("How many levels to generate? "));
  const minDifficulty = parseIntStrict(await rl.question("Min difficulty (1-10)? "));
  const maxDifficulty = parseIntStrict(await rl.question("Max difficulty (1-10)? "));

  if (!Number.isInteger(amount) || amount <= 0) {
    throw new Error("Amount must be a positive integer.");
  }
  validateRange(minDifficulty, maxDifficulty);

  const startLevel = getNextLevelNumber(levelsDir);
  const endLevel = startLevel + amount - 1;

  console.log(`\nNext level number: ${startLevel}`);
  console.log(`Generating ${amount} levels (${startLevel}..${endLevel})`);
  console.log(`Difficulty span: ${minDifficulty} to ${maxDifficulty}\n`);

  for (let levelNumber = startLevel; levelNumber <= endLevel; levelNumber++) {
    let tries = 0;

    // Keep trying until THIS level number is generated and saved.
    while (true) {
      tries++;

      // IMPORTANT: re-roll difficulty on every try, otherwise you can get stuck retrying a cursed value forever.
      const difficulty = randomInt(minDifficulty, maxDifficulty);
      const filename = `Level ${levelNumber} ${difficulty}.json`;
      const filepath = path.join(levelsDir, filename);

      const startedAt = Date.now();
      let lastProgress = null;

      console.log(`Generating level ${levelNumber} (difficulty ${difficulty}) [try ${tries}]...`);

      const heartbeat = setInterval(() => {
        const sec = Math.floor((Date.now() - startedAt) / 1000);
        const tail = lastProgress
          ? ` | attempt ${lastProgress.attempt}/${lastProgress.maxAttempts} bestΔ=${lastProgress.bestDelta} bestDiff=${lastProgress.bestEstimatedDifficulty}`
          : "";
        console.log(`  ...still working (${sec}s)${tail}`);
      }, 2000);

      try {
        const result = await generateLevel({
          requestedDifficulty: difficulty,
          levelNumber,
          onProgress: (p) => {
            lastProgress = p;
            console.log(
              `  progress: ${p.attempt}/${p.maxAttempts} | bestΔ=${p.bestDelta} | bestDiff=${p.bestEstimatedDifficulty} | unique=${p.bestUnique} | forcedFirst=${p.bestForcedFirst ?? ""}`
            );
          },
        });

        clearInterval(heartbeat);

        const levelObj = result?.level ?? result;
        const json = JSON.stringify(levelObj, null, 2);
        fs.writeFileSync(filepath, json, "utf8");

        const size = fs.statSync(filepath).size;
        const totalSec = Math.floor((Date.now() - startedAt) / 1000);

        console.log(`Saved: ${filename} (${size} bytes) ✅ (${totalSec}s)\n`);
        break;
      } catch (err) {
        clearInterval(heartbeat);

        console.log(`  Failed try ${tries}: ${err?.message ?? String(err)}`);
        const backoff = Math.min(8000, 500 + tries * 300);
        console.log(`  Retrying after ${backoff}ms...\n`);
        await sleep(backoff);
      }
    }
  }

  console.log("All levels generated. 🎯");
} catch (err) {
  console.error("\nFatal error:", err?.message ?? err);
} finally {
  rl.close();
}
