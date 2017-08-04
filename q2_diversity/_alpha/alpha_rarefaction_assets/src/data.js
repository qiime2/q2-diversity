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

let data = null;
class Data {
  constructor() {
    if (!data) {
      data = this;
    }
    this.series = {};
    return data;
  }
  appendSeries(name, series) {
    this.series[name] = series;
    this.series[name].dots = false;
    this.series[name].line = false;
    console.log('append series: ', name, series);
  }
  toggle(name, dots, line) {
    this.series[name].dots = dots;
    this.series[name].line = line;
  }
  getSeries(name) {
    return this.series[name];
  }
}
export const curData = new Data();
