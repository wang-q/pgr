//! Export FAS block variations (substitutions/indels) to an Excel workbook.

use rust_xlsxwriter::*;
use std::cmp::max;
use std::collections::BTreeMap;

use crate::libs::alignment::{
    get_indels, get_subs, polarize_indels, polarize_subs, Indel, Substitution,
};
use crate::libs::fmt::fas::{next_fas_block, FasBlock};

/// Export variations from FAS blocks to an Excel xlsx file.
#[allow(clippy::too_many_arguments)]
pub fn export_to_xlsx(
    infiles: &[String],
    outfile: &str,
    wrap: u16,
    is_indel: bool,
    is_outgroup: bool,
    no_single: bool,
    no_complex: bool,
    min_freq: Option<f64>,
    max_freq: Option<f64>,
) -> anyhow::Result<()> {
    let mut workbook = Workbook::new();
    let worksheet = workbook.add_worksheet();

    let format_of: BTreeMap<String, Format> = create_formats();

    let mut opt = Opt {
        sec_cursor: 1,
        col_cursor: 1,
        sec_height: 0,
        max_name_len: 1,
        wrap,
        color_loop: 15,
        seq_count: 0,
        is_outgroup,
    };

    for infile in infiles {
        let mut reader = crate::reader(infile)?;

        while let Ok(block) = next_fas_block(&mut reader) {
            let mut seqs: Vec<&[u8]> = vec![];
            for entry in &block.entries {
                seqs.push(entry.seq().as_ref());
            }

            let vars = get_vars(
                &seqs,
                is_outgroup,
                is_indel,
                no_single,
                no_complex,
                min_freq,
                max_freq,
            )?;

            opt.seq_count = seqs.len() as u32;
            opt.sec_height = opt.seq_count + 2;
            opt.col_cursor = 1;

            paint_name(worksheet, &format_of.clone(), &mut opt, &block)?;

            if opt.is_outgroup {
                opt.seq_count -= 1;
            }

            for (_, var) in vars {
                match var {
                    Variation::Substitution(sub) => {
                        paint_sub(worksheet, &format_of.clone(), &mut opt, &sub).unwrap()
                    }
                    Variation::Indel(indel) => {
                        paint_indel(worksheet, format_of.clone(), &mut opt, &indel)?
                    }
                }

                opt.col_cursor += 1;
                if opt.col_cursor > opt.wrap {
                    opt.col_cursor = 1;
                    opt.sec_cursor += 1;
                }
            }

            opt.sec_cursor += 1;
        }
    }

    worksheet.set_column_width(0, opt.max_name_len as f64)?;
    for i in 1..=(opt.wrap + 3) {
        worksheet.set_column_width(i, 1.6)?;
    }

    workbook.save(outfile)?;
    Ok(())
}

#[derive(Debug)]
enum Variation {
    Substitution(Substitution),
    Indel(Indel),
}

#[derive(Debug)]
struct Opt {
    sec_cursor: u32,
    col_cursor: u16,
    sec_height: u32,
    max_name_len: usize,
    wrap: u16,
    color_loop: u32,
    seq_count: u32,
    is_outgroup: bool,
}

fn paint_name(
    worksheet: &mut Worksheet,
    format_of: &BTreeMap<String, Format>,
    opt: &mut Opt,
    block: &FasBlock,
) -> anyhow::Result<()> {
    for i in 1..=block.entries.len() {
        let pos_row = opt.sec_height * (opt.sec_cursor - 1);

        let rg = block.entries[i - 1].range().to_string();
        worksheet.write_with_format(
            pos_row + i as u32,
            0,
            rg.clone(),
            format_of.clone().get("name").unwrap(),
        )?;

        opt.max_name_len = max(rg.len(), opt.max_name_len);
    }
    Ok(())
}

fn paint_indel(
    worksheet: &mut Worksheet,
    format_of: BTreeMap<String, Format>,
    opt: &mut Opt,
    indel: &Indel,
) -> anyhow::Result<()> {
    let mut pos_row = opt.sec_height * (opt.sec_cursor - 1);

    let col_taken = indel.length.min(3) as u16;

    if opt.col_cursor + col_taken > opt.wrap {
        opt.col_cursor = 1;
        opt.sec_cursor += 1;
        pos_row = opt.sec_height * (opt.sec_cursor - 1);
    }

    let indel_string = format!("{}{}", indel.itype, indel.length);
    let format = {
        let bg_idx = if indel.occurred == "unknown" {
            "unknown".to_string()
        } else {
            let idx = u32::from_str_radix(&indel.occurred, 2)? % opt.color_loop;
            idx.to_string()
        };
        let format_key = format!("indel_{}", bg_idx);
        format_of.get(&format_key).unwrap()
    };

    for i in 1..=opt.seq_count {
        let mut flag_draw = false;
        if indel.occurred == "unknown" {
            flag_draw = true;
        } else {
            let occ = indel.occurred.chars().nth(i as usize - 1).unwrap();
            if occ == '1' {
                flag_draw = true;
            }
        }

        if !flag_draw {
            continue;
        }

        if col_taken == 1 {
            worksheet.write_with_format(
                pos_row,
                opt.col_cursor,
                indel.start,
                format_of.get("pos").unwrap(),
            )?;
            worksheet.write_with_format(pos_row + i, opt.col_cursor, &indel_string, format)?;
        } else if col_taken == 2 {
            worksheet.write_with_format(
                pos_row,
                opt.col_cursor,
                indel.start,
                format_of.get("pos").unwrap(),
            )?;
            worksheet.write_with_format(
                pos_row,
                opt.col_cursor + 1,
                indel.end,
                format_of.get("pos").unwrap(),
            )?;
            worksheet.merge_range(
                pos_row + i,
                opt.col_cursor,
                pos_row + i,
                opt.col_cursor + 1,
                &indel_string,
                format,
            )?;
        } else {
            worksheet.write_with_format(
                pos_row,
                opt.col_cursor,
                indel.start,
                format_of.get("pos").unwrap(),
            )?;
            worksheet.write_with_format(
                pos_row,
                opt.col_cursor + 1,
                "|",
                format_of.get("pos").unwrap(),
            )?;
            worksheet.write_with_format(
                pos_row,
                opt.col_cursor + 2,
                indel.end,
                format_of.get("pos").unwrap(),
            )?;
            worksheet.merge_range(
                pos_row + i,
                opt.col_cursor,
                pos_row + i,
                opt.col_cursor + 2,
                &indel_string,
                format,
            )?;
        }
    }
    Ok(())
}

fn paint_sub(
    worksheet: &mut Worksheet,
    format_of: &BTreeMap<String, Format>,
    opt: &mut Opt,
    sub: &Substitution,
) -> anyhow::Result<()> {
    let pos_row = opt.sec_height * (opt.sec_cursor - 1);

    worksheet.write_with_format(
        pos_row,
        opt.col_cursor,
        sub.pos,
        format_of.get("pos").unwrap(),
    )?;

    for i in 1..=opt.seq_count {
        let base = sub.bases.chars().nth(i as usize - 1).unwrap();
        let occurred = if sub.pattern == "unknown" {
            '0'
        } else {
            sub.pattern.chars().nth(i as usize - 1).unwrap()
        };

        let base_color = if occurred == '1' {
            let bg_idx = u32::from_str_radix(&sub.pattern, 2)? % opt.color_loop;
            format!("sub_{}_{}", base, bg_idx)
        } else {
            format!("sub_{}_unknown", base)
        };
        let format = format_of.get(&base_color).unwrap();
        worksheet.write_with_format(pos_row + i, opt.col_cursor, base.to_string(), format)?;
    }

    if opt.is_outgroup {
        let base_color = format!("sub_{}_unknown", sub.obase);
        let format = format_of.get(&base_color).unwrap();
        worksheet.write_with_format(
            pos_row + opt.seq_count + 1,
            opt.col_cursor,
            sub.obase.clone(),
            format,
        )?;
    }
    Ok(())
}

fn get_vars(
    seqs: &[&[u8]],
    is_outgroup: bool,
    is_indel: bool,
    no_single: bool,
    no_complex: bool,
    min_freq: Option<f64>,
    max_freq: Option<f64>,
) -> anyhow::Result<BTreeMap<i32, Variation>> {
    let mut vars = BTreeMap::new();

    let mut seq_count = seqs.len();
    let out_seq = if is_outgroup {
        seq_count -= 1;
        Some(seqs[seq_count])
    } else {
        None
    };

    let subs = if is_outgroup {
        let mut unpolarized = get_subs(&seqs[..seq_count])?;
        polarize_subs(&mut unpolarized, out_seq.unwrap());
        unpolarized
    } else {
        get_subs(seqs)?
    };

    for sub in subs {
        if no_single && sub.freq <= 1 {
            continue;
        }
        if no_complex && sub.freq == -1 {
            continue;
        }
        if let Some(min) = min_freq {
            if (sub.freq as f64) / (seq_count as f64) < min {
                continue;
            }
        }
        if let Some(max) = max_freq {
            if (sub.freq as f64) / (seq_count as f64) > max {
                continue;
            }
        }

        vars.insert(sub.pos, Variation::Substitution(sub));
    }

    if is_indel {
        let indels = if is_outgroup {
            let mut unpolarized = get_indels(&seqs[..seq_count])?;
            polarize_indels(&mut unpolarized, out_seq.unwrap())?;
            unpolarized
        } else {
            get_indels(seqs)?
        };

        for indel in indels {
            if no_single && indel.freq <= 1 {
                continue;
            }
            if no_complex && indel.freq == -1 {
                continue;
            }
            if let Some(min) = min_freq {
                if (indel.freq as f64) / (seq_count as f64) < min {
                    continue;
                }
            }
            if let Some(max) = max_freq {
                if (indel.freq as f64) / (seq_count as f64) > max {
                    continue;
                }
            }

            vars.insert(indel.start, Variation::Indel(indel));
        }
    }

    Ok(vars)
}

fn create_formats() -> BTreeMap<String, Format> {
    let mut format_of: BTreeMap<String, Format> = BTreeMap::new();

    format_of.insert(
        "name".to_string(),
        Format::new().set_font_name("Courier New").set_font_size(10),
    );

    format_of.insert(
        "pos".to_string(),
        Format::new()
            .set_font_name("Courier New")
            .set_font_size(8)
            .set_align(FormatAlign::VerticalCenter)
            .set_align(FormatAlign::Center)
            .set_rotation(90),
    );

    let bg_colors: Vec<u32> = vec![
        0xC0C0C0, 0xFFFF99, 0xCCFFCC, 0xCCFFFF, 0x99CCFF, 0xCC99FF, 0xFFCC99, 0x9999FF, 0x33CCCC,
        0xFFCC00, 0xFF99CC, 0xFF9900, 0xFFFFCC, 0xFF8080, 0xCCCCFF,
    ];

    let sub_fc_of: BTreeMap<String, u32> = BTreeMap::from([
        ("A".to_string(), 0x003300),
        ("C".to_string(), 0x000080),
        ("G".to_string(), 0x660066),
        ("T".to_string(), 0x800000),
        ("N".to_string(), 0x000000),
        ("N".to_string(), 0x000000),
        ("-".to_string(), 0x000000),
    ]);

    for fc in sub_fc_of.keys() {
        format_of.insert(
            format!("sub_{}_{}", fc, "unknown"),
            Format::new()
                .set_font_name("Courier New")
                .set_font_size(10)
                .set_align(FormatAlign::VerticalCenter)
                .set_align(FormatAlign::Center)
                .set_font_color(*sub_fc_of.get(fc).unwrap())
                .set_background_color(Color::White),
        );

        for i in 0..bg_colors.len() {
            let key = format!("sub_{}_{}", fc, i);
            format_of.insert(
                key,
                Format::new()
                    .set_font_name("Courier New")
                    .set_font_size(10)
                    .set_align(FormatAlign::VerticalCenter)
                    .set_align(FormatAlign::Center)
                    .set_font_color(*sub_fc_of.get(fc).unwrap())
                    .set_background_color(*bg_colors.get(i).unwrap()),
            );
        }
    }

    format_of.insert(
        "sub_-".to_string(),
        Format::new()
            .set_font_name("Courier New")
            .set_font_size(10)
            .set_align(FormatAlign::VerticalCenter)
            .set_align(FormatAlign::Center),
    );

    for i in 0..bg_colors.len() {
        let key = format!("indel_{}", i);
        format_of.insert(
            key,
            Format::new()
                .set_font_name("Courier New")
                .set_font_size(10)
                .set_bold()
                .set_align(FormatAlign::VerticalCenter)
                .set_align(FormatAlign::Center)
                .set_background_color(*bg_colors.get(i).unwrap()),
        );
    }
    format_of.insert(
        format!("indel_{}", "unknown"),
        Format::new()
            .set_font_name("Courier New")
            .set_font_size(10)
            .set_bold()
            .set_align(FormatAlign::VerticalCenter)
            .set_align(FormatAlign::Center)
            .set_background_color(Color::White),
    );

    format_of
}
