/* global d */

import render from './render'; // warning, stats

export default function setupData(data, metric) {
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

function updateData(metric, category, svg, href) {
  href.attr('href', `${metric}.csv`);
  let data = d[metric];
  if (category) {
    data = d[metric][category];
  }
  const preppedData = setupData(data, metric);
  render(svg, preppedData, category);
}

let curState = null;
class State {
  constructor() {
    if (!curState) {
      curState = this;
    }
    this.category = '';
    this.metric = '';
    this.svg = null;
    this.href = null;
    return curState;
  }
  initialize(metric, category, row, svg) {
    // CONTROLS
    const downloadDiv = row.append('div')
      .attr('class', 'col-lg-2 form-group downloadCSV');
    this.href = downloadDiv.append('a')
        .attr('href', '')
        .text('Download CSV');
    this.svg = svg;
    this.metric = metric;
    this.category = category;
    updateData(metric, category, this.svg, this.href);
  }
  setCategory(c) {
    this.category = c;
    updateData(this.metric, this.category, this.svg, this.href);
  }
  setMetric(m) {
    this.metric = m;
    updateData(this.metric, this.category, this.svg, this.href);
  }
  getCategory() {
    return this.category;
  }
  getMetric() {
    return this.metric;
  }
}
export const state = new State();
