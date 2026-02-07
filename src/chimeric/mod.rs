// Chimeric alignment detection module
//
// Detects split reads that span distant genomic locations:
// - Inter-chromosomal fusions (e.g., BCR-ABL)
// - Intra-chromosomal strand breaks
// - Circular RNAs (back-splices)
//
// Detection strategy:
// - Tier 1: Soft-clip analysis (post-alignment)
// - Tier 2: Multi-cluster chimeric stitching (during clustering)

mod detect;
mod output;
mod score;
mod segment;

pub use detect::ChimericDetector;
pub use output::ChimericJunctionWriter;
pub use segment::{ChimericAlignment, ChimericSegment};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_module_exports() {
        // Ensure all public types are accessible
        let _ = std::mem::size_of::<ChimericAlignment>();
        let _ = std::mem::size_of::<ChimericSegment>();
    }
}
