from flask import Flask, request, jsonify
from cfg_converter import parse_grammar_from_lines, CFGtoChomsky

app = Flask(__name__)


@app.after_request
def add_cors_headers(response):
    # Allow requests from local files or other origins during local dev
    response.headers['Access-Control-Allow-Origin'] = '*'
    response.headers['Access-Control-Allow-Methods'] = 'GET, POST, OPTIONS'
    response.headers['Access-Control-Allow-Headers'] = 'Content-Type'
    return response

def serialize_productions(productions):
    """Serialize productions without renaming symbols.
    - Deduplicate RHS strings per LHS
    - Represent empty RHS as 'ε'
    - Keep original variable names (e.g., T0, T1, A1) to match deterministic CNF output
    """
    out = {}
    for lhs in sorted(productions.keys()):
        seen = set()
        lst = []
        for rhs in productions.get(lhs, []):
            if not rhs:
                rep = 'ε'
            else:
                rep = ''.join(rhs)
            if rep not in seen:
                seen.add(rep)
                lst.append(rep)
        out[lhs] = lst
    return out

@app.route('/convert', methods=['POST'])
def convert():
    data = request.get_json() or {}
    grammar_text = data.get('grammar', '')
    if not grammar_text:
        return jsonify({'error':'no grammar provided'}), 400

    lines = grammar_text.splitlines()
    cfg = parse_grammar_from_lines(lines)
    converter = CFGtoChomsky(cfg)

    steps = {}
    # Step 1: add new start
    converter._add_new_start(verbose=False)
    steps['after_add_start'] = serialize_productions(converter.cfg.productions)

    # Step 2: eliminate epsilon
    converter._eliminate_epsilon(verbose=False)
    steps['after_epsilon'] = serialize_productions(converter.cfg.productions)

    # Step 3: eliminate unit productions
    converter._eliminate_unit(verbose=False)
    steps['after_unit'] = serialize_productions(converter.cfg.productions)

    # Step 4: convert to cnf
    converter._convert_to_cnf(verbose=False)
    steps['final'] = serialize_productions(converter.cfg.productions)

    return jsonify({
        'steps': steps,
        'start_symbol': converter.cfg.start_symbol,
    })

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)
