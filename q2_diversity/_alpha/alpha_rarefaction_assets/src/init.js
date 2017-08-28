/* global categories */
/* global metrics */
import { select } from 'd3';
import state from './state';
import { addMetricPicker, addCategoryPicker } from './toolbar';

export default function init() {
  const metric = metrics[0];
  const category = categories[0];
  const body = select('#main');
  const controlsRow = body.select('.controls');
  const plotSvg = body.select('.plotSvg');
  const legendBox = body.select('.legendBoxSvg');
  const legendTitle = body.select('.legendTitle');
  state.initialize(metric, category, controlsRow, plotSvg, legendBox, legendTitle);
  addMetricPicker(controlsRow, metrics, metric);
  if (categories.length > 0) {
    addCategoryPicker(controlsRow, categories, category);
  }
}
