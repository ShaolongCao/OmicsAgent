"""Export a checked public source tree, excluding local/private artifacts."""
import argparse
import shutil
from pathlib import Path

from release_check import scan_tree

ROOT_FILES = ['README.md','pyproject.toml','uv.lock','.gitignore','.gitattributes']
DOCS = ['HIERARCHY_TUTORIAL.md','methods.md','skills.md','validation-v03.md','RELEASE_NOTES.md',
        'SCIENTIST_REVIEW.md', 'validation-review.md']


def export(root: Path, destination: Path) -> None:
    if destination.exists():
        raise ValueError('Release destination must not exist')
    selected = [root/name for name in ROOT_FILES] + [root/'docs'/name for name in DOCS]
    selected += [p for p in (root/'docs'/'cases').rglob('*') if p.is_file()]
    for folder in ['src','tests','examples','scripts','.github']:
        selected += [p for p in (root/folder).rglob('*') if p.is_file()
                     and '__pycache__' not in p.parts and p.suffix != '.pyc']
    for source in selected:
        if source.is_symlink():
            raise ValueError('Symlink in release source')
        target = destination / source.relative_to(root)
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(source,target)
    findings = scan_tree(destination)
    if findings:
        raise ValueError(f'Release checks failed (paths/rules only): {findings}')
    print(f'Exported {len(selected)} checked files to {destination}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,default=Path('.'))
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    export(args.root.resolve(),args.out.resolve())
