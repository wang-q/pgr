use std::io::Write;

use super::NamedMatrix;

/// PHYLIP matrix output format.
pub enum MatrixFormat {
    Full,
    Lower,
    Strict,
}

impl MatrixFormat {
    pub fn from_mode(s: &str) -> anyhow::Result<Self> {
        match s {
            "full" => Ok(Self::Full),
            "lower" => Ok(Self::Lower),
            "strict" => Ok(Self::Strict),
            _ => anyhow::bail!("unsupported output format"),
        }
    }
}

/// Write a NamedMatrix in the specified PHYLIP format.
pub fn write_phylip_matrix<W: Write>(
    m: &NamedMatrix,
    fmt: MatrixFormat,
    writer: &mut W,
) -> anyhow::Result<()> {
    let names = m.get_names();
    let size = m.size();

    writer.write_fmt(format_args!("{:>4}\n", size))?;

    for (i, name) in names.iter().enumerate().take(size) {
        match fmt {
            MatrixFormat::Full => {
                writer.write_fmt(format_args!("{}", name))?;
                for j in 0..size {
                    writer.write_fmt(format_args!("\t{}", m.get(i, j)))?;
                }
            }
            MatrixFormat::Lower => {
                writer.write_fmt(format_args!("{}", name))?;
                for j in 0..i {
                    writer.write_fmt(format_args!("\t{}", m.get(i, j)))?;
                }
            }
            MatrixFormat::Strict => {
                writer.write_fmt(format_args!(
                    "{:<10}",
                    name.chars().take(10).collect::<String>()
                ))?;
                for j in 0..size {
                    writer.write_fmt(format_args!(" {:.6}", m.get(i, j)))?;
                }
            }
        }
        writer.write_fmt(format_args!("\n"))?;
    }

    Ok(())
}

/// Write a submatrix restricted to `names`. Returns the list of names not found in `m`.
pub fn write_subset<W: Write>(
    m: &NamedMatrix,
    names: &[String],
    writer: &mut W,
) -> anyhow::Result<Vec<String>> {
    let all_names = m.get_names();
    let mut indices = Vec::new();
    let mut missing = Vec::new();

    for name in names {
        match m.get_index(name) {
            Some(idx) => indices.push(idx),
            None => missing.push(name.clone()),
        }
    }

    writer.write_fmt(format_args!("{}\n", indices.len()))?;

    for &i in &indices {
        writer.write_fmt(format_args!("{}", all_names[i]))?;
        for &j in &indices {
            writer.write_fmt(format_args!("\t{}", m.get(i, j)))?;
        }
        writer.write_fmt(format_args!("\n"))?;
    }

    Ok(missing)
}

/// Extract paired values from the upper triangle (excluding diagonal) of two matrices,
/// restricted to sequence names common to both. Returns `(common_names, values1, values2)`.
pub fn extract_common_upper_triangle(
    m1: &NamedMatrix,
    m2: &NamedMatrix,
) -> anyhow::Result<(Vec<String>, Vec<f32>, Vec<f32>)> {
    let names1 = m1.get_names();
    let names2 = m2.get_names();
    let common_names: Vec<String> = names1
        .iter()
        .filter(|name| names2.contains(name))
        .map(|s| s.to_string())
        .collect();

    if common_names.is_empty() {
        anyhow::bail!("No common sequence names found between matrices");
    }

    let mut values1 = Vec::with_capacity(common_names.len() * (common_names.len() - 1) / 2);
    let mut values2 = Vec::with_capacity(common_names.len() * (common_names.len() - 1) / 2);

    for i in 0..common_names.len() {
        for j in 0..i {
            if let (Some(v1), Some(v2)) = (
                m1.get_by_name(&common_names[i], &common_names[j]),
                m2.get_by_name(&common_names[i], &common_names[j]),
            ) {
                values1.push(v1);
                values2.push(v2);
            }
        }
    }

    Ok((common_names, values1, values2))
}
