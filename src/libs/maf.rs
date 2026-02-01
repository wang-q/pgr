use std::io::{self, Write};

#[derive(Debug, Clone, Default)]
pub struct MafComp {
    pub src: String,
    pub start: usize,
    pub size: usize,
    pub strand: char,
    pub src_size: usize,
    pub text: String,
}

#[derive(Debug, Clone, Default)]
pub struct MafAli {
    pub score: Option<f64>,
    pub components: Vec<MafComp>,
}

pub struct MafWriter<W: Write> {
    writer: W,
}

impl<W: Write> MafWriter<W> {
    pub fn new(writer: W) -> Self {
        Self { writer }
    }

    pub fn write_header(&mut self, program: &str) -> io::Result<()> {
        writeln!(self.writer, "##maf version=1 scoring={}", program)
    }

    pub fn write_ali(&mut self, ali: &MafAli) -> io::Result<()> {
        writeln!(self.writer, "a score={:.1}", ali.score.unwrap_or(0.0))?;
        for comp in &ali.components {
            writeln!(
                self.writer,
                "s {:<20} {:10} {:10} {} {:10} {}",
                comp.src, comp.start, comp.size, comp.strand, comp.src_size, comp.text
            )?;
        }
        writeln!(self.writer)?;
        Ok(())
    }
}
