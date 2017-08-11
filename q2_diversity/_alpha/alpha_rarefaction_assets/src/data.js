import state from './state';

export function setupData(data, metric) {
  const [xAxisLabel, yAxisLabel] = ['Sequencing Depth', metric];
  let minX = Infinity;
  let maxX = 0;
  let minY = Infinity;
  let maxY = 0;
  let minSubY = Infinity;
  let maxSubY = 0;
  const depthIndex = data.columns.indexOf('depth');
  const medianIndex = data.columns.indexOf('50%');
  const countIndex = data.columns.indexOf('count');
  data.data.forEach((d) => {
    const x = d[depthIndex];
    const y = d[medianIndex];
    if (x < minX) minX = x;
    if (x > maxX) maxX = x;
    if (y < minY) minY = y;
    if (y > maxY) maxY = y;
    const count = d[countIndex];
    console.log(count);
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
      .selectAll(`[class="circle ${name}"]`)
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
