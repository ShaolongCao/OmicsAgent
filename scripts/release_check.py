"""Check public artifacts without disclosing credential matches."""
import argparse
import json
import re
from pathlib import Path

RULES = {
    'cjk': re.compile(r'[\u3400-\u9fff\uf900-\ufaff]'),
    'absolute_home': re.compile(r'/(?:Users|home)/[^/\s]+/'),
    'credential': re.compile(r'sk-ant-[A-Za-z0-9_-]{20,}|sk-[A-Za-z0-9_-]{24,}|'
                             r'gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,}|'
                             r'-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----|'
                             r'https?://[^/\s:@]+:[^/\s@]+@'),
}
FORBIDDEN = {'private_source','data','outputs','.venv','.cache','.superpowers','__pycache__'}


def scan_tree(root: Path) -> list[dict[str, str]]:
    findings = []
    for path in sorted(root.rglob('*')):
        relative = path.relative_to(root)
        if '.git' in relative.parts:
            continue
        if path.is_symlink():
            findings.append({'file':str(relative),'rule':'symlink'})
            continue
        if not path.is_file():
            continue
        if (set(relative.parts) & FORBIDDEN or path.name.startswith('.env')
                or path.suffix in {'.h5ad','.pem','.key','.p12'}):
            findings.append({'file':str(relative),'rule':'forbidden_path'})
            continue
        try:
            content = path.read_text()
        except UnicodeError:
            findings.append({'file':str(relative),'rule':'unreviewed_binary'})
            continue
        for name, pattern in RULES.items():
            if pattern.search(content):
                findings.append({'file':str(relative),'rule':name})
    return findings


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('root',type=Path)
    args = parser.parse_args()
    findings = scan_tree(args.root)
    print(json.dumps({'findings':findings,'passed':not findings},indent=2))
    raise SystemExit(bool(findings))
