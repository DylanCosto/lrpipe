#!/usr/bin/env python3
"""Copy the license files of every crate compiled into the release binaries into a bundle.

usage: tools/third_party_licenses.py OUT_DIR [cargo build flags...]

For each package outside the workspace that the release build links (normal dependencies of
the workspace's binaries, for the host platform and these flags), copies its LICENSE*, LICENCE*,
COPYING*, NOTICE* and COPYRIGHT* files to OUT_DIR/NAME-VERSION/, and writes OUT_DIR/INDEX.tsv
with each package's name, version, license expression and source. A package that ships no
license file is listed with its expression and `no license file`, and the script fails if one
has no license expression either, so a release never ships a dependency without its terms.
"""
import json
import pathlib
import shutil
import subprocess
import sys

PREFIXES = ('license', 'licence', 'copying', 'notice', 'copyright')


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    out = pathlib.Path(sys.argv[1])
    flags = sys.argv[2:]
    host = next(line.split()[1] for line in subprocess.check_output(['rustc', '-vV'], text=True).splitlines()
                if line.startswith('host:'))
    metadata = json.loads(subprocess.check_output(
        ['cargo', 'metadata', '--format-version', '1', '--locked', '--filter-platform', host, *flags]))
    packages = {package['id']: package for package in metadata['packages']}
    members = set(metadata['workspace_members'])
    nodes = {node['id']: node for node in metadata['resolve']['nodes']}
    # Everything reachable from the workspace through normal (not dev or build) dependencies.
    linked, pending = set(), list(members)
    while pending:
        node = nodes[pending.pop()]
        for dependency in node['deps']:
            if any(kind['kind'] is None for kind in dependency['dep_kinds']) and dependency['pkg'] not in linked:
                linked.add(dependency['pkg'])
                pending.append(dependency['pkg'])
    out.mkdir(parents=True, exist_ok=True)
    rows, missing = [], []
    for package_id in sorted(linked - members, key=lambda p: (packages[p]['name'], packages[p]['version'])):
        package = packages[package_id]
        name, version = package['name'], package['version']
        root = pathlib.Path(package['manifest_path']).parent
        files = sorted(f for f in root.iterdir() if f.is_file() and f.name.lower().startswith(PREFIXES))
        if package.get('license_file'):
            extra = root / package['license_file']
            if extra.is_file() and extra not in files:
                files.append(extra)
        expression = package.get('license') or ''
        if not files and not expression:
            missing.append(f'{name} {version}')
            continue
        target = out / f'{name}-{version}'
        target.mkdir(exist_ok=True)
        for f in files:
            shutil.copy(f, target / f.name)
        source = package.get('source') or 'path'
        rows.append(f"{name}\t{version}\t{expression or '-'}\t{source}\t{'' if files else 'no license file'}")
    if missing:
        sys.exit('no license file or license expression: ' + ', '.join(missing))
    (out / 'INDEX.tsv').write_text('name\tversion\tlicense\tsource\tnote\n' + '\n'.join(rows) + '\n')
    print(f'{len(rows)} third-party packages: licenses in {out}')


if __name__ == '__main__':
    main()
