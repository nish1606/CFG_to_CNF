# CFG to CNF Converter

A small toolkit to convert a Context-Free Grammar (CFG) into Chomsky Normal Form (CNF). It includes:

- A Python CLI (`cfg_converter.py`) that performs stepwise conversion: add start (if needed), ε-elimination, unit-elimination, terminal factoring in mixed rules, and binarization.
- A simple static webpage (`interactive_cnf.html`) that parses an input grammar and shows the conversion steps and the final grouped CNF.
- An optional Flask server endpoint (if you add `app.py`) that returns conversion steps as JSON for the UI. (The page can run fully client-side without a server.)

## Features

- Strict CNF: final productions are either A → a or A → BC; ε only allowed for the start symbol if the language includes ε.
- Selective terminal replacement: terminals are only factored out when they appear in longer productions.
- Deterministic helper naming for readability (e.g., X0, X1…).
- Final grouped display: one line per nonterminal with `|` separators and readable spacing.

## Folder layout

```
.
├─ cfg_converter.py        # Python CLI tool for CFG→CNF
├─ interactive_cnf.html    # Client-side CNF visualizer (open in a browser)
└─ sample_grammar.txt      # Example grammars (input format)
```

> Note: If you created a Flask backend (`app.py`), it is optional; the HTML works without it.

## Grammar format

- One production per line: `A → aB | bC | ε`
- You can use either `→` or `->` for the arrow.
- Separate alternatives with `|`. A `/` is also accepted by the tools and treated as `|`.
- Nonterminals: uppercase symbols like `S`, `A`, `B`.
- Terminals: lowercase letters or digits (e.g., `a`, `b`, `0`, `1`).
- `ε` denotes the empty string.

Examples:

```
S → AB | a
A → aA | ε
B → bB | b
```

```
S -> a S a | b S b | a | b | ε
```

## Using the Python CLI

Run the converter against a grammar file or standard input.

- From a file:

```powershell
python d:\CFG_to_CNF\cfg_converter.py -f d:\CFG_to_CNF\sample_grammar.txt
```

- From stdin (paste your grammar, then Ctrl+Z and Enter on Windows):

```powershell
python d:\CFG_to_CNF\cfg_converter.py
```

Useful flags:
- `--no-verbose` to reduce intermediate logging.
- `--trace` to print a detailed, educational step-by-step trace (for small grammars).
- `--canonical` to print a textbook-style, deterministic CNF trace.

## Using the web page (no backend required)

1. Start a simple static server in this folder (or open the file directly):

```powershell
python -m http.server 8000
```

2. Open the page in your browser:

```
http://127.0.0.1:8000/interactive_cnf.html
```

3. Paste your grammar into the Input area and click "Convert to CNF". The page shows:
   - Original grammar
   - ε-elimination
   - Unit-elimination
   - Final CNF (grouped)

The page also verifies CNF constraints and allows copying/downloading the final grammar.

### Optional: convert via server

If you add a Flask backend (example `app.py`) with a `/convert` endpoint, the page can also call the server. Click "Convert via Server" to use it. This is optional; all conversion can run client-side.

## Examples

`sample_grammar.txt` contains commented examples. Uncomment one group at a time when testing. Example:

```
S → ASB
A → aAS | a | ε 
B → SbS | A | bb
```

## Troubleshooting

- Empty output:
  - Ensure your input lines aren’t commented out with `#`.
  - Use `→` or `->`, and separate alternatives with `|` (or `/`).
- Epsilon handling:
  - ε is removed from all rules except possibly the start symbol if nullable.
- Unit productions:
  - A → B rules are eliminated via closure so A inherits B’s non-unit rules.
- Strict CNF:
  - Single terminals stay as A → a; terminals inside longer RHS are factored to helper variables first.
- Browser cache:
  - Hard refresh (Ctrl+F5) if the page appears unchanged after edits.

## License

This project is provided as-is for educational use. You’re free to adapt it for coursework or personal projects.
