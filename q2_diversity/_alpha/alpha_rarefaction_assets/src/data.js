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

  data.data.data.forEach((d) => {
    console.log(d);
    const x = d.depth;
    const vals = d.slice(2);
    const yMax = max(vals);
    const yMin = min(vals);
    minX = min(x, minX);
    maxX = max(x, maxX);
    minY = min(yMin, minY);
    maxY = max(yMax, maxY);
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
