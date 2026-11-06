import json

from test_review import fixture_run


def test_portable_bundle_preserves_measurements_and_records_redactions(tmp_path):
    from omicsagent.review_bundle import export_review_bundle
    (tmp_path/'source').mkdir()
    source = fixture_run(tmp_path/'source')
    private_root = 'C:' + '/' + 'Users/private'
    (source/'config.json').write_text(json.dumps({'input': private_root+'/query.h5ad'}))
    (source/'manifest.json').write_text(json.dumps({'run_id': 'example', 'artifacts': {}}))
    (source/'trace.jsonl').write_text(json.dumps({'path': private_root+'/query.h5ad'})+'\n')
    destination = tmp_path/'portable'
    export_review_bundle(source, destination)
    assert (destination/'hierarchy'/'decisions.csv').read_bytes() == (source/'hierarchy'/'decisions.csv').read_bytes()
    assert private_root not in (destination/'config.json').read_text()
    assert private_root not in (destination/'trace.jsonl').read_text()
    manifest = json.loads((destination/'manifest.json').read_text())
    assert manifest['status'] == 'derived_review_bundle'
    assert manifest['redacted_files'] == ['config.json', 'trace.jsonl']
    assert (destination/'report.html').is_file()
    assert private_root in (source/'config.json').read_text()
