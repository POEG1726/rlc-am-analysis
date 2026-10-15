"""Check publishable files, archive hashes and common privacy leaks without hardware."""
import ast
import csv
import hashlib
import json
from pathlib import Path
import re
import struct
import subprocess

ROOT = Path(__file__).resolve().parents[1]
SKIP = {'.git','.venv','venv','__pycache__','outputs','local_runs','.pytest_cache','.idea','.vscode'}
TEXT = {'.py','.md','.json','.csv','.txt','.yml','.yaml','.toml','.html','.svg'}
PATTERNS = {
    'device resource or unique identifier': re.compile(r'USB\d*::|(?:SDS|SDG)[0-9][A-Z0-9]{7,}'),
    'machine absolute path': re.compile(r'(?<![A-Za-z0-9])[A-Za-z]:[\\/]|/(?:home|Users)/[^\s/]+'),
    'email address': re.compile(r'[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}'),
    'assigned credential': re.compile(r'''(?i)(?:api[_-]?key|password|access[_-]?token|client[_-]?secret)\s*[=:]\s*["'][^"']+["']'''),
}


def publishable_files():
    if (ROOT/'.git').exists():
        result=subprocess.run(['git','ls-files','--cached','--others','--exclude-standard','-z'],
                              cwd=ROOT,check=True,capture_output=True)
        return sorted({ROOT / s.decode('utf-8') for s in result.stdout.split(b'\0') if s})
    return sorted(p for p in ROOT.rglob('*') if p.is_file() and not SKIP.intersection(p.relative_to(ROOT).parts)
                  and not p.name.endswith('.local.json'))


def main():
    errors=[];files=publishable_files();pictures=0;python_files=0
    for p in files:
        rel=p.relative_to(ROOT).as_posix()
        if SKIP.intersection(p.relative_to(ROOT).parts) or p.suffix in {'.npz','.log','.pyc'} or p.name=='scope_test.csv' or p.name.endswith('.local.json') or p.name.startswith('.env'):
            errors.append(f'{rel}: private/generated file is publishable')
        if p.stat().st_size > 5*1024**2:errors.append(f'{rel}: unexpected large file')
        if p.suffix in TEXT or p.name in {'.gitignore','.gitattributes'}:
            text=p.read_text(encoding='utf-8-sig')
            for title,pattern in PATTERNS.items():
                if pattern.search(text):errors.append(f'{rel}: possible {title}')
            if p.suffix=='.py':
                ast.parse(text,filename=rel);python_files+=1
            if p.suffix=='.json':json.loads(text)
            if p.suffix=='.csv':
                with p.open(encoding='utf-8',newline='') as f:
                    rows=csv.reader(f);header=next(rows)
                    if any(len(r)!=len(header) for r in rows):errors.append(f'{rel}: malformed CSV row')
        if p.suffix=='.png':
            pictures+=1;data=p.read_bytes();offset=8
            while offset+12<=len(data):
                size=struct.unpack('>I',data[offset:offset+4])[0];kind=data[offset+4:offset+8]
                value=data[offset+8:offset+8+size];offset+=size+12
                if kind in {b'iTXt',b'zTXt',b'eXIf'} or (kind==b'tEXt' and not value.startswith(b'Software\0Matplotlib')):
                    errors.append(f'{rel}: unexpected embedded image metadata')
    archive=json.loads((ROOT/'data/manifest.json').read_text(encoding='utf-8'))
    for item in archive['files']:
        p=ROOT/item['file']
        if hashlib.sha256(p.read_bytes()).hexdigest()!=item['sha256'] or p.stat().st_size!=item['bytes']:
            errors.append(f'{item["file"]}: archive hash or size mismatch')
        if item['rows'] is not None:
            with p.open(encoding='utf-8',newline='') as f:
                if sum(1 for _ in csv.DictReader(f))!=item['rows']:errors.append(f'{item["file"]}: archive row count mismatch')
    example=json.loads((ROOT/'config.example.json').read_text())
    if example['scope_address'] or example['generator_address']:
        errors.append('Example configuration must have empty instrument addresses')
    if pictures!=43:errors.append(f'Expected 43 preserved PNGs; found {pictures}')
    if errors:raise SystemExit('\n'.join(errors))
    print(f'PASS: {len(files)} publishable files, {python_files} Python files, {pictures} preserved PNGs; '
          f'{len(archive["files"])} archive hashes and CSV row counts verified. No common private identifiers or absolute paths found.')


if __name__=='__main__':main()
