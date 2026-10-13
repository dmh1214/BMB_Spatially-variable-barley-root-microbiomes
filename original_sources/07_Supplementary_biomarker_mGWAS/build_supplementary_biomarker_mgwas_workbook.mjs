import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const scriptDir = path.dirname(fileURLToPath(import.meta.url));
const baseDir = path.dirname(scriptDir);
const outputDir = path.join(scriptDir, "outputs");
const csvPath = path.join(scriptDir, "supplementary_biomarker_mgwas_summary.csv");
const outPath = path.join(outputDir, "Supplementary_Table_Biomarker_mGWAS_summary.xlsx");

const csvText = await fs.readFile(csvPath, "utf8");
const workbook = await Workbook.fromCSV(csvText, { sheetName: "Biomarker Summary" });
const summary = workbook.worksheets.getItem("Biomarker Summary");
const notes = workbook.worksheets.add("README");

summary.showGridLines = false;
summary.freezePanes.freezeRows(1);
summary.freezePanes.freezeColumns(4);

const used = summary.getUsedRange();
const header = summary.getRange("A1:Z1");
header.format.fill = { color: "#DCE6F2" };
header.format.font = { bold: true, color: "#1F1F1F" };
header.format.borders = { preset: "all", style: "thin", color: "#B8C7D9" };
header.format.rowHeight = 34;
header.format.wrapText = true;

used.format.font = { name: "Aptos", size: 10 };
used.format.borders = { preset: "inside", style: "thin", color: "#E6EAF0" };
summary.getRange("A1:Z69").format.borders = {
  top: { style: "thin", color: "#B8C7D9" },
  bottom: { style: "thin", color: "#B8C7D9" },
  left: { style: "thin", color: "#B8C7D9" },
  right: { style: "thin", color: "#B8C7D9" },
  insideHorizontal: { style: "thin", color: "#E6EAF0" },
  insideVertical: { style: "thin", color: "#E6EAF0" },
};

summary.getRange("A:Z").format.wrapText = false;
summary.getRange("E:J").setNumberFormat("0.0000");
summary.getRange("I:J").setNumberFormat("0.00E+00");
summary.getRange("O:R").setNumberFormat("0");
summary.getRange("T:T").setNumberFormat("0.00E+00");
summary.getRange("Y:Y").setNumberFormat("0.00E+00");

summary.getRange("A:A").format.columnWidthPx = 110;
summary.getRange("B:B").format.columnWidthPx = 78;
summary.getRange("C:C").format.columnWidthPx = 82;
summary.getRange("D:D").format.columnWidthPx = 180;
summary.getRange("E:J").format.columnWidthPx = 90;
summary.getRange("K:L").format.columnWidthPx = 96;
summary.getRange("M:N").format.columnWidthPx = 88;
summary.getRange("O:R").format.columnWidthPx = 88;
summary.getRange("S:T").format.columnWidthPx = 90;
summary.getRange("U:U").format.columnWidthPx = 96;
summary.getRange("V:V").format.columnWidthPx = 88;
summary.getRange("W:W").format.columnWidthPx = 140;
summary.getRange("X:X").format.columnWidthPx = 88;
summary.getRange("Y:Y").format.columnWidthPx = 110;
summary.getRange("Z:Z").format.columnWidthPx = 260;

summary.getRange("M2:M69").conditionalFormats.add("containsText", {
  text: "Yes",
  format: { fill: { color: "#E8F1FB" }, font: { color: "#0B5394" } },
});
summary.getRange("N2:N69").conditionalFormats.add("containsText", {
  text: "Yes",
  format: { fill: { color: "#D9EAD3" }, font: { bold: true, color: "#274E13" } },
});
summary.getRange("N2:N69").conditionalFormats.add("containsText", {
  text: "No",
  format: { fill: { color: "#FDE9E7" }, font: { color: "#7F0000" } },
});
summary.getRange("L2:L69").conditionalFormats.add("containsText", {
  text: "Yes",
  format: { fill: { color: "#D9EAD3" }, font: { color: "#274E13" } },
});
summary.getRange("L2:L69").conditionalFormats.add("containsText", {
  text: "Inconsistent",
  format: { fill: { color: "#FCE5CD" }, font: { color: "#7F6000", bold: true } },
});
summary.getRange("W2:W69").conditionalFormats.add("containsText", {
  text: "Yes",
  format: { fill: { color: "#E8F1FB" }, font: { color: "#0B5394" } },
});

const notesRows = [
  ["Supplementary biomarker-mGWAS summary", ""],
  ["Description", "Each row represents one biomass-associated biomarker identified by feature selection in one field trial."],
  ["SHAP_sign / LME_sign / Direction_match", "SHAP-derived direction, LME-derived direction, and whether the two directions are consistent for a given biomarker."],
  ["mGWAS_tested", "Yes if this biomarker was included in the chromosome-level mGWAS summary files."],
  ["mGWAS_sig", "Yes if the biomarker had at least one significant SNP in the chromosome-level mGWAS summary tables."],
  ["In_network_summary", "Yes if the biomarker also appeared in Network/mgwas_sum.xlsx, which stores SNP-level hits used in the network summary."],
  ["Taxon", "Reported genus-level or nearest available taxonomic annotation from the feature-selection output tables."],
  ["Min_p_chr / Min_chr_p", "Chromosome with the smallest minimum p-value across all tested chromosomes and that minimum p-value, retained even when mGWAS was not significant."],
  ["Network_lead_SNP / Network_SNP_list", "Lead SNP and SNP list only when available in the network-level mGWAS summary workbook."],
  ["Source files", "Bac_Aut/feature_selection/*.csv; Fun_Aut/feature_selection/*.csv; Bac_Aut/mGWAS/*_ASV_chr_summary.tsv; Fun_Aut/mGWAS/*_OTU_chr_summary.tsv; Network/mgwas_sum.xlsx"],
];

notes.getRange(`A1:B${notesRows.length}`).values = notesRows;
notes.showGridLines = false;
notes.getRange("A1:B1").merge();
notes.getRange("A1").format.fill = { color: "#DCE6F2" };
notes.getRange("A1").format.font = { bold: true, size: 13 };
notes.getRange(`A2:B${notesRows.length}`).format.borders = {
  preset: "all",
  style: "thin",
  color: "#D9D9D9",
};
notes.getRange("A2:A9").format.font = { bold: true };
notes.getRange("A:A").format.columnWidthPx = 185;
notes.getRange("B:B").format.columnWidthPx = 560;
notes.getRange(`A1:B${notesRows.length}`).format.wrapText = true;
notes.getRange(`A1:B${notesRows.length}`).format.autofitRows();

await fs.mkdir(outputDir, { recursive: true });
const xlsx = await SpreadsheetFile.exportXlsx(workbook);
await xlsx.save(outPath);

const preview = await workbook.render({
  sheetName: "Biomarker Summary",
  range: "A1:N20",
  scale: 2,
  format: "png",
});
const previewBytes = new Uint8Array(await preview.arrayBuffer());
await fs.writeFile(path.join(outputDir, "Supplementary_Table_Biomarker_mGWAS_summary_preview.png"), previewBytes);

const inspect = await workbook.inspect({
  kind: "table",
  range: "A1:Z10",
  include: "values",
  tableMaxRows: 10,
  tableMaxCols: 26,
});
console.log(inspect.ndjson);
console.log(`Saved workbook to ${outPath}`);
