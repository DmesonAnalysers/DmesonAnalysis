## Downloading train outputs
In an environment connected to GRID (e.g. O2Physics), you can download the train outputs using the following command, after setting the correct configuration in `train_downloader.sh`:
```bash
bash train_downloader.sh
```

For MC, `parquet` files with different species and origins can be extracted using the `MC_SELECTIONS` variable. For example, for a $D_s$ meson analysis, you may use
```bash
MC_SELECTIONS=(
    "PromptDs:abs(fFlagMcMatchRec) == 4 and fOriginMcRec == 1 and fFlagMcDecayChanRec == 1"
    "PromptDplus:abs(fFlagMcMatchRec) == 4 and fOriginMcRec == 1 and fFlagMcDecayChanRec == 3"
    "Prompt:abs(fFlagMcMatchRec) == 4 and fOriginMcRec == 1 and (fFlagMcDecayChanRec == 3 or fFlagMcDecayChanRec == 1)"
    "NonPromptDs:abs(fFlagMcMatchRec) == 4 and fOriginMcRec == 2 and fFlagMcDecayChanRec == 1"
    "NonPromptDplus:abs(fFlagMcMatchRec) == 4 and fOriginMcRec == 2 and fFlagMcDecayChanRec == 3"
    "NonPrompt:abs(fFlagMcMatchRec) == 4 and fOriginMcRec == 2 and (fFlagMcDecayChanRec == 3 or fFlagMcDecayChanRec == 1)"
    "DplusBkg:abs(fFlagMcMatchRec) == 1"
    )
```
to separate prompt and non-prompt $D_s$ and $D^+$ mesons, prompt and non-prompt $D^+$ mesons.

