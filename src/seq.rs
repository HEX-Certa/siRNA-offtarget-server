//! Nucleotide helpers: normalize, reverse-complement, 2-bit packing, Hamming.

pub const MIN_OLIGO: usize = 19;
pub const MAX_OLIGO: usize = 23;
pub const SEED_START: usize = 1;
pub const SEED_LEN: usize = 7;
pub const MAX_MISMATCH: u8 = 3;
/// siDirect near-match window: oligo positions 2–20 (1-based) → 19-mer.
pub const OLIGO_2_20_LEN: usize = 19;

#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub enum SeqError {
    #[error("sequence must be {MIN_OLIGO}–{MAX_OLIGO} nt, got {0}")]
    BadLength(usize),
    #[error("sequence contains invalid base {0:?}")]
    BadBase(char),
    #[error("sequence is empty")]
    Empty,
}

/// Trim, uppercase, U→T. Does not validate alphabet.
pub fn normalize(seq: &str) -> String {
    seq.chars()
        .filter(|c| !c.is_whitespace())
        .map(|c| match c {
            'u' | 'U' => 'T',
            c => c.to_ascii_uppercase(),
        })
        .collect()
}

pub fn validate_oligo(seq: &str) -> Result<(), SeqError> {
    if seq.is_empty() {
        return Err(SeqError::Empty);
    }
    if !(MIN_OLIGO..=MAX_OLIGO).contains(&seq.len()) {
        return Err(SeqError::BadLength(seq.len()));
    }
    for c in seq.chars() {
        if !matches!(c, 'A' | 'C' | 'G' | 'T') {
            return Err(SeqError::BadBase(c));
        }
    }
    Ok(())
}

pub fn revcomp(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'N' => 'N',
            other => other,
        })
        .collect()
}

/// Guide positions 2–8 (1-based) → 0-based [1, 8).
pub fn seed_of_guide(guide: &str) -> Result<&str, SeqError> {
    if guide.len() < SEED_START + SEED_LEN {
        return Err(SeqError::BadLength(guide.len()));
    }
    Ok(&guide[SEED_START..SEED_START + SEED_LEN])
}

/// Oligo positions 2–20 (1-based) as a 19-mer. A 19-nt input is used as-is.
pub fn oligo_2_20(seq: &str) -> Result<&str, SeqError> {
    if seq.len() < OLIGO_2_20_LEN {
        return Err(SeqError::BadLength(seq.len()));
    }
    if seq.len() == OLIGO_2_20_LEN {
        return Ok(seq);
    }
    // 1-based [2, 20] → 0-based [1, 20)
    Ok(&seq[1..20])
}

/// Encode A/C/G/T → 0..=3. Anything else is N (returns 0, `is_n = true`).
#[inline]
pub fn encode_base(b: u8) -> (u8, bool) {
    match b {
        b'A' | b'a' => (0, false),
        b'C' | b'c' => (1, false),
        b'G' | b'g' => (2, false),
        b'T' | b't' | b'U' | b'u' => (3, false),
        _ => (0, true),
    }
}

#[inline]
pub fn decode_base(bits: u8) -> u8 {
    match bits & 3 {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        _ => b'T',
    }
}

/// Pack up to 32 bases into a u64; low 2 bits are the first base.
pub fn pack_oligo(seq: &str) -> u64 {
    let mut v = 0u64;
    for (i, b) in seq.bytes().take(32).enumerate() {
        let (bits, _) = encode_base(b);
        v |= (bits as u64) << (2 * i);
    }
    v
}

/// Hamming distance of two 2-bit packed oligos (first `len` bases).
#[inline]
pub fn hamming_packed(a: u64, b: u64, len: usize) -> u32 {
    debug_assert!(len <= 32);
    let xor = a ^ b;
    let diffs = (xor | (xor >> 1)) & 0x5555_5555_5555_5555;
    let mask = if len >= 32 {
        u64::MAX
    } else {
        (1u64 << (2 * len)) - 1
    };
    (diffs & mask).count_ones()
}

pub fn hamming_str(a: &str, b: &str) -> u32 {
    a.bytes().zip(b.bytes()).filter(|(x, y)| x != y).count() as u32
}

/// Extract `len` bases starting at `pos` from a 2-bit packed stream (32 bases / u64).
#[inline]
pub fn extract_packed(packed: &[u64], pos: usize, len: usize) -> u64 {
    debug_assert!(len <= 32);
    let bit = pos * 2;
    let word = bit / 64;
    let shift = bit % 64;
    let bits_needed = len * 2;
    let mut v = packed[word] >> shift;
    if shift + bits_needed > 64 {
        v |= packed[word + 1] << (64 - shift);
    }
    if bits_needed >= 64 {
        v
    } else {
        v & ((1u64 << bits_needed) - 1)
    }
}

/// True if any N is set in `start..start+len` of a 1-bit-per-base mask (64 pos / u64).
#[inline]
pub fn window_has_n(nmask: &[u64], start: usize, len: usize) -> bool {
    let end = start + len;
    let mut pos = start;
    while pos < end {
        let w = pos / 64;
        let bit = pos % 64;
        let take = (end - pos).min(64 - bit);
        let mask = if take == 64 {
            u64::MAX
        } else {
            ((1u64 << take) - 1) << bit
        };
        if nmask[w] & mask != 0 {
            return true;
        }
        pos += take;
    }
    false
}

pub fn decode_window(packed: &[u64], pos: usize, len: usize) -> String {
    let bits = extract_packed(packed, pos, len);
    let mut s = String::with_capacity(len);
    for i in 0..len {
        s.push(decode_base(((bits >> (2 * i)) & 3) as u8) as char);
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalize_u_and_whitespace() {
        assert_eq!(normalize(" auC gU "), "ATCGT");
    }

    #[test]
    fn revcomp_guide_example() {
        let guide = "ATAAACTCCAGGCCTATGAGG";
        assert_eq!(revcomp(guide), "CCTCATAGGCCTGGAGTTTAT");
    }

    #[test]
    fn seed_slice() {
        let guide = "ATAAACTCCAGGCCTATGAGG";
        assert_eq!(seed_of_guide(guide).unwrap(), "TAAACTC");
        assert_eq!(revcomp("TAAACTC"), "GAGTTTA");
    }

    #[test]
    fn oligo_2_20_21mer_and_19mer() {
        let guide21 = "ATAAACTCCAGGCCTATGAGG";
        assert_eq!(oligo_2_20(guide21).unwrap(), "TAAACTCCAGGCCTATGAG");
        assert_eq!(oligo_2_20(guide21).unwrap().len(), 19);
        let guide19 = "TAAACTCCAGGCCTATGAG";
        assert_eq!(oligo_2_20(guide19).unwrap(), guide19);
        assert!(matches!(
            oligo_2_20("ATAA"),
            Err(SeqError::BadLength(4))
        ));
    }

    #[test]
    fn hamming_packed_matches_str() {
        let a = "ATAAACTCCAGGCCTATGAGG";
        let b = "ATAAACTCCAGGCCTATGAGA";
        assert_eq!(hamming_packed(pack_oligo(a), pack_oligo(b), 21), 1);
        assert_eq!(hamming_str(a, b), 1);
        assert_eq!(hamming_packed(pack_oligo(a), pack_oligo(a), 21), 0);
    }

    #[test]
    fn extract_roundtrip() {
        let seq = "ACGTACGTAAATAAACTCCAGGCCTATGAGGTTTT";
        let mut packed = vec![0u64; 4];
        for (i, b) in seq.bytes().enumerate() {
            let (bits, _) = encode_base(b);
            packed[i / 32] |= (bits as u64) << (2 * (i % 32));
        }
        assert_eq!(decode_window(&packed, 10, 21), &seq[10..31]);
    }

    #[test]
    fn validate_rejects_bad() {
        assert!(validate_oligo("ATAAACTCCAGGCCTATGA").is_ok());
        assert!(matches!(
            validate_oligo("ATAA"),
            Err(SeqError::BadLength(4))
        ));
        assert!(matches!(
            validate_oligo("ATAAACTCCAGGCCTATGANN"),
            Err(SeqError::BadBase('N'))
        ));
    }
}
