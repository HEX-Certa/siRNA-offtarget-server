use crate::types::OfftargetCounts;

/// T6b specificity label from non-target counts.
pub fn specificity_label(off: &OfftargetCounts) -> &'static str {
    if off.perfect > 0 || off.mismatch_1 > 0 {
        "低"
    } else if off.mismatch_2 > 0 || off.mismatch_3 > 0 || off.seed_transcripts > 0 {
        "中"
    } else {
        "高"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(p: u64, m1: u64, m2: u64, m3: u64, seed: u64) -> OfftargetCounts {
        OfftargetCounts {
            perfect: p,
            mismatch_1: m1,
            mismatch_2: m2,
            mismatch_3: m3,
            seed_transcripts: seed,
            seed_genes: seed,
        }
    }

    #[test]
    fn labels() {
        assert_eq!(specificity_label(&c(0, 0, 0, 0, 0)), "高");
        assert_eq!(specificity_label(&c(0, 0, 2, 0, 0)), "中");
        assert_eq!(specificity_label(&c(0, 0, 0, 1, 0)), "中");
        assert_eq!(specificity_label(&c(0, 0, 0, 0, 5)), "中");
        assert_eq!(specificity_label(&c(1, 0, 0, 0, 0)), "低");
        assert_eq!(specificity_label(&c(0, 1, 0, 0, 0)), "低");
        assert_eq!(specificity_label(&c(0, 1, 3, 0, 9)), "低");
    }
}
