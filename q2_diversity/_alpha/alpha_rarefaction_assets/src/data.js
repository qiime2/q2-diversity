import state from './state';

export function setupData(data, metric) {
  const [xAxisLabel, yAxisLabel] = ['Sequencing Depth', metric];
  let minX = Infinity;
  let maxX = 0;
  let minY = Infinity;
  let maxY = 0;
  const depthIndex = data.columns.indexOf('depth');
  const minIndex = data.columns.indexOf('min');
  const maxIndex = data.columns.indexOf('max');
  data.data.forEach((d) => {
    const x = d[depthIndex];
    if (x < minX) minX = x;
    if (x > maxX) maxX = x;
    if (d[minIndex] < minY) minY = d[minIndex];
    if (d[maxIndex] > maxY) maxY = d[maxIndex];
  });

  return {
    data,
    xAxisLabel,
    yAxisLabel,
    minX,
    maxX,
    minY,
    maxY,
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
