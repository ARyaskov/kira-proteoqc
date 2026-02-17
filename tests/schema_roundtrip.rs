use kira_proteoqc::schema::v1::{Mode, ProteoQcV1};

#[test]
fn schema_roundtrip_v1() {
    let report = ProteoQcV1::empty("0.0.0-test", Mode::Cell, false, true);
    let json = serde_json::to_string(&report).unwrap();
    let decoded: ProteoQcV1 = serde_json::from_str(&json).unwrap();
    assert_eq!(decoded.tool, "kira-proteoqc");
    assert_eq!(decoded.schema_version, "v1");
    assert!(matches!(decoded.input_meta.mode, Mode::Cell));
}
