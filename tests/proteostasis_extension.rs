use kira_proteoqc::metrics::proteostasis_extension::panels::{
    AGGREGATION_PANEL, CHAPERONE_PANEL, ERAD_PANEL, PROTEASOME_PANEL, PROTEO_EXTENSION_PANEL_V1,
    UPR_PANEL,
};

#[test]
fn proteostasis_panel_contract_is_stable() {
    assert_eq!(PROTEO_EXTENSION_PANEL_V1, "proteo_extension_panel_v1");
    assert_eq!(CHAPERONE_PANEL.len(), 8);
    assert_eq!(PROTEASOME_PANEL.len(), 10);
    assert_eq!(UPR_PANEL.len(), 6);
    assert_eq!(ERAD_PANEL.len(), 4);
    assert_eq!(AGGREGATION_PANEL.len(), 3);
}
