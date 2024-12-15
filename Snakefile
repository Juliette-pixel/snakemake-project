configfile: "config.yaml"

rule all:
    input:
        config["outputs"]["umap_plot"].format(output_dir=config["output_dir"]),
        config["outputs"]["cluster_markers"].format(output_dir=config["output_dir"])

rule setup_and_initialization:
    input:
        data_dir=config["data_dir"]
    output:
        config["outputs"]["initial"].format(output_dir=config["output_dir"])
    params:
        script=config["scripts"]["setup"]
    shell:
        """
        Rscript {params.script} {input.data_dir} {output}
        """

rule qc_and_preprocessing:
    input:
        config["outputs"]["initial"].format(output_dir=config["output_dir"])
    output:
        config["outputs"]["preprocessed"].format(output_dir=config["output_dir"])
    params:
        script=config["scripts"]["qc"]
    shell:
        """
        Rscript {params.script} {input} {output}
        """

rule normalization_and_feature_selection:
    input:
        config["outputs"]["preprocessed"].format(output_dir=config["output_dir"])
    output:
        config["outputs"]["normalized"].format(output_dir=config["output_dir"])
    params:
        script=config["scripts"]["normalization"]
    shell:
        """
        Rscript {params.script} {input} {output}
        """

rule scaling_and_pca:
    input:
        config["outputs"]["normalized"].format(output_dir=config["output_dir"])
    output:
        config["outputs"]["pca"].format(output_dir=config["output_dir"])
    params:
        script=config["scripts"]["scaling"]
    shell:
        """
        Rscript {params.script} {input} {output}
        """

rule clustering_and_umap:
    input:
        config["outputs"]["pca"].format(output_dir=config["output_dir"])
    output:
        config["outputs"]["umap"].format(output_dir=config["output_dir"])
    params:
        script=config["scripts"]["clustering"]
    shell:
        """
        Rscript {params.script} {input} {output}
        """

rule differential_markers_and_visualization:
    input:
        config["outputs"]["umap"].format(output_dir=config["output_dir"])
    output:
        markers=config["outputs"]["cluster_markers"].format(output_dir=config["output_dir"]),
        plot=config["outputs"]["umap_plot"].format(output_dir=config["output_dir"])
    params:
        script=config["scripts"]["markers"]
    shell:
        """
        Rscript {params.script} {input} {output.markers} {output.plot}
        """
