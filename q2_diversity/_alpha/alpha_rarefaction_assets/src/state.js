/* global d */
import { setupData } from './data';
import render from './render';

function updateData(metric, column, svg, href, legend, legendTitle) {
  href.attr('href', `${metric}.csv`);
  let data = d[metric];
  if (column) {
    data = d[metric][column];
  }
  const preppedData = setupData(data, metric);
  render(svg, preppedData, column, legend, legendTitle);
}

class State {
  constructor() {
    this.column = '';
    this.metric = '';
    this.svg = null;
    this.href = null;
    this.legend = null;
    this.legendTitle = null;
  }
  initialize(metric, column, row, svg, legend, legendTitle) {
    this.href = row.select('.downloadCSV a');
    this.svg = svg;
    this.metric = metric;
    this.column = column;
    this.legend = legend;
    this.legendTitle = legendTitle;
    updateData(metric, column, this.svg, this.href, this.legend, this.legendTitle);
  }
  setColumn(c) {
    this.column = c;
    updateData(this.metric, this.column, this.svg, this.href, this.legend, this.legendTitle);
  }
  setMetric(m) {
    this.metric = m;
    updateData(this.metric, this.column, this.svg, this.href, this.legend, this.legendTitle);
  }
  getColumn() {
    return this.column;
  }
  getMetric() {
    return this.metric;
  }
  getSvg() {
    return this.svg;
  }
}
const state = new State();
export { state as default };
