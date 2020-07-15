import state from './state';

export function setupData(data, metric) {
  const [xAxisLabel, yAxisLabel] = ['Sequencing Depth', metric];
  let minX = Infinity;
  let maxX = 0;
  let minY = Infinity;
  let maxY = 0;
  let minSubY = Infinity;
  let maxSubY = 0;
  const depthIndex = data.columns.indexOf('_alpha_rarefaction_depth_column_');
  const lowerWhiskerIndex = data.columns.indexOf('9%');
  const upperWhiskerIndex = data.columns.indexOf('91%');
  const countIndex = data.columns.indexOf('count');
  data.data.forEach((d) => {
    const x = d[depthIndex];
    if (x < minX) minX = x;
    if (x > maxX) maxX = x;
    const yMin = d[lowerWhiskerIndex];
    const yMax = d[upperWhiskerIndex];
    if (yMin < minY) minY = yMin;
    if (yMax > maxY) maxY = yMax;
    const count = d[countIndex];
    if (count > maxSubY) maxSubY = count;
    if (count < minSubY) minSubY = count;
  });

  return {
    data,
    xAxisLabel,
    yAxisLabel,
    minX,
    maxX,
    minY,
    maxY,
    minSubY,
    maxSubY,
  };
}


const curData = {};
export function appendSeries(name, series, curColor) {
  curData[name] = series;
  curData[name].dotsOpacity = 1;
  curData[name].lineOpacity = 1;
  curData[name].dots = curColor;
  curData[name].line = curColor;
}
export function toggle(name, dots, line) {
  if (dots !== null) {
    curData[name].dots = dots;
    curData[name].dotsOpacity = dots === 'white' ? 0 : 1;
    // update chart
    state.getSvg()
      .selectAll(`[class="symbol ${name}"]`)
      .attr('opacity', curData[name].dotsOpacity);
  }
  if (line !== null) {
    curData[name].line = line;
    curData[name].lineOpacity = line === 'white' ? 0 : 1;
    // update chart
    state.getSvg()
      .selectAll(`[class="line ${name}"]`)
      .attr('opacity', curData[name].lineOpacity);
  }
}
export { curData };
