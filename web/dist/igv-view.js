import igv from "/vendor/igv-3.8.9.esm.min.js";

const status = document.querySelector("#igv-status");
const runId = new URLSearchParams(location.search).get("run");

async function loadBrowser() {
  const response = await fetch(`/api/runs/${encodeURIComponent(runId)}/igv-view`);
  if (!response.ok) throw new Error("IGV 轨道数据暂不可用");
  const {loci, tracks} = await response.json();
  const signal = {type:"wig", height:55, autoscaleGroup:"CPM", min:0};
  const peaks = {type:"annotation", height:34, displayMode:"COLLAPSED", indexed:false};
  const browser = await igv.createBrowser(document.querySelector("#igv-browser"), {
    loadDefaultGenomes:false,
    reference:{id:"hg38", name:"Human GRCh38/hg38", format:"chromsizes", url:"/api/igv/hg38.chrom.sizes"},
    locus:loci[0] || "chr1:1000000-1010000",
    showIdeogram:false,
    showSVGButton:false,
    tracks:tracks.map(track => ({...(track.format === "bigwig" ? signal : peaks), ...track})),
  });
  const regions = document.querySelector("#igv-loci");
  for (const locus of loci) {
    const button = document.createElement("button");
    button.type = "button";
    button.textContent = locus;
    button.addEventListener("click", () => browser.search(locus));
    regions.append(button);
  }
  regions.hidden = !loci.length;
  status.hidden = true;
}

loadBrowser().catch(error => {
  status.textContent = `轨道加载失败。可使用 IGV Desktop 会话查看。 ${error.message}`;
});
