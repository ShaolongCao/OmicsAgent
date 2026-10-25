import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location('release_check',Path(__file__).parents[1]/'scripts'/'release_check.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_detects_secrets_without_returning_values(tmp_path):
    secret = 'ghp_' + 'A'*36
    (tmp_path/'bad.txt').write_text(secret)
    findings = module.scan_tree(tmp_path)
    assert findings == [{'file':'bad.txt','rule':'credential'}]
    assert secret not in repr(findings)


def test_detects_cjk_and_private_paths(tmp_path):
    (tmp_path/'bad.md').write_text('\u4e2d\u6587 ' + '/' + 'Users/example/private')
    assert {x['rule'] for x in module.scan_tree(tmp_path)} == {'cjk','absolute_home'}


def test_rejects_environment_files_and_symlinks(tmp_path):
    (tmp_path/'.env').write_text('KEY=placeholder')
    (tmp_path/'link').symlink_to(tmp_path/'.env')
    assert {x['rule'] for x in module.scan_tree(tmp_path)} == {'forbidden_path','symlink'}
