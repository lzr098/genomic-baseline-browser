"""Minimal Chinese overlays for variant report PDFs (titles + glossary definitions)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

PdfLocale = str  # "en" | "zh"

NOTO_CJK_REGULAR = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
NOTO_CJK_BOLD = "/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc"
NOTO_CJK_SC_INDEX = 2


@dataclass(frozen=True)
class PdfTitles:
    section_variant_overview: str
    section_annotation_interpretation: str
    section_frequency_transcripts: str
    section_sample_comparison: str
    section_sample_baseline_comparison: str
    section_population_frequency: str
    section_data_sources: str
    section_glossary: str
    subsection_annotation_scores: str
    subsection_clinvar: str
    subsection_interpretation: str
    subsection_sample_finding: str
    subsection_baseline_verdict: str


TITLES_EN = PdfTitles(
    section_variant_overview="VARIANT OVERVIEW",
    section_annotation_interpretation="ANNOTATION & INTERPRETATION",
    section_frequency_transcripts="FREQUENCY & TRANSCRIPTS",
    section_sample_comparison="SAMPLE & COMPARISON",
    section_sample_baseline_comparison="SAMPLE & BASELINE COMPARISON",
    section_population_frequency="POPULATION FREQUENCY",
    section_data_sources="DATA SOURCES",
    section_glossary="GLOSSARY — ANNOTATION & SCORES",
    subsection_annotation_scores="ANNOTATION & SCORES",
    subsection_clinvar="CLINVAR",
    subsection_interpretation="INTERPRETATION SUMMARY",
    subsection_sample_finding="SAMPLE FINDING",
    subsection_baseline_verdict="BASELINE VERDICT",
)

TITLES_ZH = PdfTitles(
    section_variant_overview="变体概述",
    section_annotation_interpretation="注释&解释",
    section_frequency_transcripts="频率与转录本",
    section_sample_comparison="样本与基线对比",
    section_sample_baseline_comparison="样本与基线对比",
    section_population_frequency="人群频率",
    section_data_sources="数据来源",
    section_glossary="术语表 — 注释与评分",
    subsection_annotation_scores="注释与评分",
    subsection_clinvar="CLINVAR",
    subsection_interpretation="解读摘要",
    subsection_sample_finding="样本检测结果",
    subsection_baseline_verdict="基线判定",
)

GLOSSARY_DEFINITIONS_ZH: dict[str, str] = {
    "Gene": "变异位点对应的 HGNC 基因符号；若有全称则在括号中给出。",
    "Function": "来自基因数据库（如 MyGene.info / NCBI Gene）的简要功能介绍。",
    "Transcript": "用于主要 HGVS 命名的规范 Ensembl 转录本（优先 MANE Select）。",
    "Consequence": "描述对转录本影响的 Sequence Ontology 术语（如 intron_variant）。",
    "Impact": "VEP 预测的功能影响严重程度：HIGH、MODERATE、LOW 或 MODIFIER。",
    "Exon": "受影响外显子编号（当前/总数）；若不在外显子内则为 n/a。",
    "Intron": "内含子变异对应的内含子编号（当前/总数）。",
    "HGVS c.": "相对于所列转录本的 HGVS 编码 DNA 描述。",
    "HGVS p.": "HGVS 蛋白改变；若无蛋白水平改变则为 n/a。",
    "Constr. tx": "当与规范转录本不同时，用于 gnomAD v4 基因约束指标的转录本。",
    "pLI": "功能缺失不耐受概率（gnomAD）；接近 1 表示基因对 LoF 耐受较差。",
    "LOEUF": "LoF 观察/期望上界分数（gnomAD）；越低表示 LoF 约束越强。",
    "mis_z": "错义约束 Z 分数（gnomAD）；越高表示错义变异越稀缺。",
    "CADD Phred": "整合多种注释的 CADD PHRED 分数；越高提示有害性可能越大。",
    "REVEL max": "REVEL 错义致病性最高分（0–1）；越高越支持致病性错义效应。",
    "SpliceAI max": "SpliceAI 最大 delta 分数（0–1）；≥0.5 常视为剪接影响值得关注。",
    "SIFT max": "受影响转录本中 SIFT 最高分；越低通常表示替代越不被耐受。",
    "PolyPhen max": "PolyPhen-2 预测有害错义替代的最高分。",
    "phyloP": "phyloP 核苷酸保守性最高分；正值表示进化保守。",
    "Sample call": "个体样本 VCF 在该位点检出的基因型与杂合状态。",
    "Genotype": "样本 VCF 的 GT 字段（如 0/1 表示一条 ref、一条 alt）。",
    "Zygosity": "由基因型推导的杂合状态：heterozygous、homozygous_alt 等。",
    "Depth": "样本 VCF 中该位点的测序深度（FORMAT DP）。",
    "Allele depth": "支持 ref 与 alt 的 reads 数（FORMAT AD：ref / alt）。",
    "Quality": "样本 VCF 记录的质量分（QUAL 或样本 GQ）。",
    "FILTER": "样本 VCF 的 FILTER 状态（如 PASS 或未通过滤器标签）。",
    "Match status": "相对 gnomAD exomes + ClinVar 基线的机器可读分类（如 known_gnomad_clinvar、novel_in_sample）。",
    "Verdict": "样本等位基因与基线关系的简短可读摘要。",
    "Novel": "若在 gnomAD exomes 与 ClinVar 基线均未出现则为 Yes（novel_in_sample）。",
    "Priority": "报告分诊用的启发式分数与等级；novel 或临床相关变异通常更高。",
    "In gnomAD exomes": "该等位基因是否存在于本地 gnomAD v4.1 exomes sites VCF（WES 基线）。",
    "In gnomAD genomes": "该等位基因是否存在于本地 gnomAD v4.1 genomes sites VCF（补充信息，不参与 novel 判定）。",
    "In ClinVar": "该等位基因是否存在于本地 ClinVar GRCh38 VCF。",
}


def get_pdf_titles(locale: PdfLocale) -> PdfTitles:
    return TITLES_ZH if locale == "zh" else TITLES_EN


def glossary_definition(term: str, definition_en: str, locale: PdfLocale) -> str:
    if locale != "zh":
        return definition_en
    return GLOSSARY_DEFINITIONS_ZH.get(term, definition_en)


def variant_report_pdf_path(json_path: str | Any, locale: PdfLocale = "en") -> Path:
    path = Path(json_path)
    name = path.name
    if name.endswith(".report.json"):
        stem = name[: -len(".report.json")]
        if locale == "en":
            return path.parent / f"{stem}.report.pdf"
        return path.parent / f"{stem}.report.{locale}.pdf"
    if locale == "en":
        return path.with_suffix(".pdf")
    return path.with_name(f"{path.stem}.{locale}.pdf")
