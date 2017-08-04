/* global d */
import setupData from './data';
import render from './render';

function updateData(metric, category, svg, href, legend, legendTitle) {
  href.attr('action', `${metric}.csv`);
  let data = d[metric];
  if (category) {
    data = d[metric][category];
  }
  const preppedData = setupData(data, metric);
  render(svg, preppedData, category, legend, legendTitle);
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
    this.legend = null;
    this.legendTitle = null;
    return curState;
  }
  initialize(metric, category, row, svg, legend, legendTitle) {
    // CONTROLS
    const downloadDiv = row.append('div')
      .attr('class', 'col-lg-2 form-group downloadCSV');
    this.href = downloadDiv.append('form')
      .attr('method', 'GET')
      .attr('action', '#');
    this.href.append('button')
      .attr('class', 'btn btn-block btn-primary btn-md form-control')
      .attr('role', 'button')
      .attr('aria-pressed', 'true')
      .text('Download CSV');
    this.svg = svg;
    this.metric = metric;
    this.category = category;
    this.legend = legend;
    this.legendTitle = legendTitle;
    updateData(metric, category, this.svg, this.href, this.legend, this.legendTitle);
  }
  setCategory(c) {
    this.category = c;
    updateData(this.metric, this.category, this.svg, this.href, this.legend, this.legendTitle);
  }
  setMetric(m) {
    this.metric = m;
    updateData(this.metric, this.category, this.svg, this.href, this.legend, this.legendTitle);
  }
  getCategory() {
    return this.category;
  }
  getMetric() {
    return this.metric;
  }
}
const state = new State();
export { state as default };
