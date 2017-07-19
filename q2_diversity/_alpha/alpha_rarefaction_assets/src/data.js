import {
  max,
  min,
} from 'd3';

export default function setupData(data, metric) {
  const [xAxisLabel, yAxisLabel] = ['Sequencing Depth', metric];

  let minX = Infinity;
  let maxX = 0;
  let minY = Infinity;
  let maxY = 0;
  console.log(data);
  const depthIndex = data.columns.indexOf('depth');
  const minIndex = data.columns.indexOf('min');
  const maxIndex = data.columns.indexOf('max');
  data.data.forEach((d) => {
    console.log('evaluating d: depth - ', d[depthIndex], ', min - ', d[minIndex], ', max - ', d[maxIndex]);
    const x = d[depthIndex];
    minX = min(x, minX);
    maxX = max(x, maxX);
    minY = min(d[minIndex], minY);
    maxY = max(d[maxIndex], maxY);
  });

  return {
    data: data.data.data,
    xAxisLabel,
    yAxisLabel,
    minX,
    maxX,
    minY,
    maxY,
  };
}

let curState = null;
class State {
  constructor() {
    if (!curState) {
      curState = this;
    }
    this.category = '';
    this.metric = '';
    return curState;
  }
  setCategory(c) {
    this.category = c;
    console.log('set to: ', this.category);
  }
  setMetric(m) {
    this.metric = m;
    console.log('set to: ', this.metric);
  }
  getCategory() {
    return this.category;
  }
  getMetric() {
    return this.metric;
  }
}
export const state = new State();
